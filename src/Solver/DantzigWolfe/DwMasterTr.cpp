/*
 * DwMasterTr.cpp
 *
 *  Created on: Nov 18, 2016
 *      Author: kibaekkim
 */

//#define DSP_DEBUG

#include "cplex.h"
/** Coin */
#include "OsiCpxSolverInterface.hpp"
/** Dsp */
#include <DantzigWolfe/DwMasterTr.h>

DSP_RTN_CODE DwMasterTr::createProblem() {
#define FREE_MEMORY \
	FREE_PTR(mat) \
	FREE_ARRAY_PTR(obj) \
	FREE_ARRAY_PTR(clbd) \
	FREE_ARRAY_PTR(cubd) \
	FREE_ARRAY_PTR(rlbd) \
	FREE_ARRAY_PTR(rubd) \
	FREE_ARRAY_PTR(starts_tr) \
	FREE_ARRAY_PTR(rows_tr) \
	FREE_ARRAY_PTR(elements_tr)

	OsiCpxSolverInterface* cpx = NULL;

	/** master problem */
	CoinPackedMatrix * mat = NULL;
	double* obj = NULL;
	double* clbd = NULL;
	double* cubd = NULL;
	double* rlbd = NULL;
	double* rubd = NULL;

	/** trust-region columns */
	int* starts_tr = NULL;
	int* rows_tr = NULL;
	double* elements_tr = NULL;

	BGN_TRY_CATCH

	/** number of initial columns */
	ncols_tr_ = 2 * (nrows_orig_ + nrows_branch_);

	/** column index of the first generated columns added */
	ncols_start_ = ncols_tr_;

	/** allocate memory */
	obj = new double [ncols_tr_];
	clbd = new double [ncols_tr_];
	cubd = new double [ncols_tr_];
	rlbd = new double [nrows_];
	rubd = new double [nrows_];
	tr_center_ = new double [nrows_orig_ + nrows_branch_];
	starts_tr = new int [ncols_tr_ + 1];
	rows_tr = new int [ncols_tr_];
	elements_tr = new double [ncols_tr_];

	/** initial trust region center */
	CoinZeroN(tr_center_, nrows_orig_ + nrows_branch_);

	/** create column-wise matrix and set number of rows */
	mat = new CoinPackedMatrix(true, 0, 0);
	mat->setDimensions(nrows_, 0);

	/** add variables related to trust region */
	tr_size_ = par_->getDblParam("DW/TR/SIZE");
	CoinFillN(obj, ncols_tr_, tr_size_);
	CoinZeroN(clbd, ncols_tr_);
	CoinFillN(cubd, ncols_tr_, COIN_DBL_MAX);
	for (int j = 0, k = 0; j < nrows_orig_ + nrows_branch_; ++j) {
		starts_tr[k] = k;
		rows_tr[k] = j;
		elements_tr[k] = 1.0;
		starts_tr[k+1] = k+1;
		rows_tr[k+1] = j;
		elements_tr[k+1] = -1.0;
		k += 2;
	}
	starts_tr[ncols_tr_] = ncols_tr_;
	mat->appendCols(ncols_tr_, starts_tr, rows_tr, elements_tr);

	DSPdebug(mat->verifyMtx(4));

	/** Set row bounds */
	CoinCopyN(org_rlbd_, nrows_orig_, rlbd);
	CoinCopyN(org_rubd_, nrows_orig_, rubd);
	CoinFillN(rlbd + nrows_orig_, nrows_conv_, 1.0);
	CoinFillN(rubd + nrows_orig_, nrows_conv_, 1.0);
	for (int i = 0; i < nrows_branch_; ++i) {
		int j = branch_row_to_col_[nrows_core_ + i];
		rlbd[nrows_core_ + i] = org_clbd_[j];
		rubd[nrows_core_ + i] = org_cubd_[j];
	}

	/** create solver */
	si_ = new OsiCpxSolverInterface();

	/** message setting */
	si_->messageHandler()->logLevel(1);

	/** load problem data */
	si_->loadProblem(*mat, clbd, cubd, obj, rlbd, rubd);

	/** write mps */
	DSPdebug(si_->writeMps("master"));

	/** NOTE: This does not go to Phase 1. */
	phase_ = 2;

	DSP_RTN_CHECK_RTN_CODE(initialOsiSolver());

	/** initial solve */
	//si_->initialSolve();

	/** TODO: what should be primal solution? in the original? or the restricted? */
	/** allocate memory for solution */
	bestprimsol_ = new double [ncols_orig_];
	primsol_ = new double [ncols_orig_];

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	/** release memory */
	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

DSP_RTN_CODE DwMasterTr::solve() {
	double stime;
	BGN_TRY_CATCH

	//DSP_RTN_CHECK_RTN_CODE(restoreCols());
	tic();

	itercnt_ = 0;
	while (1) {
		status_ = DSP_STAT_FEASIBLE;
		DSP_RTN_CHECK_RTN_CODE(solvePhase2());

		/** collect solutions */
		if (status_ == DSP_STAT_OPTIMAL ||
				status_ == DSP_STAT_FEASIBLE ||
				status_ == DSP_STAT_LIM_ITERorTIME) {
			/** recover original solution */
			CoinZeroN(primsol_, ncols_orig_);
			for (unsigned k = 0, j = ncols_tr_; k < cols_generated_.size(); ++k) {
				if (cols_generated_[k]->active_) {
					CoinPackedVector xlam = cols_generated_[k]->x_ * si_->getColSolution()[j];
					for (int i = 0; i < xlam.getNumElements(); ++i)
						primsol_[xlam.getIndices()[i]] += xlam.getElements()[i];
					j++;
				}
			}

			/** heuristics */
			DSP_RTN_CHECK_RTN_CODE(heuristics());
			break;
		} else if (status_ == DSP_STAT_DUAL_INFEASIBLE)
			break;
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DwMasterTr::solvePhase1() {
	BGN_TRY_CATCH
	DSP_RTN_CHECK_RTN_CODE(solvePhase2());
	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

DSP_RTN_CODE DwMasterTr::solvePhase2() {
#define FREE_MEMORY \
	FREE_ARRAY_PTR(price) \
	FREE_ARRAY_PTR(piA) \
	FREE_PTR(ws)

	double* price = NULL;
	double* piA = NULL;
	int nColsAdded = 0;
	double curlb = 0.0;

	CoinWarmStartBasis* ws = NULL;

	BGN_TRY_CATCH

	phase_ = 2;
	primobj_ = +COIN_DBL_MAX;
	dualobj_ = -COIN_DBL_MAX;

	/** allocate memory */
	price = new double [nrows_];
	piA = new double [ncols_orig_];

	/** initialize trust region */
	tr_size_ = par_->getDblParam("DW/TR/SIZE");

	/** initial price */
	CoinZeroN(price, nrows_);
	CoinFillN(price + nrows_orig_, nrows_conv_, COIN_DBL_MAX);

	/** generate initial columns */
	DSP_RTN_CHECK_RTN_CODE(generateCols(price, piA, curlb, nColsAdded));

	/** termination test: subproblem solution may declare infeasible. */
	if (status_ == DSP_STAT_OPTIMAL || status_ == DSP_STAT_FEASIBLE) {
		DSPdebugMessage("%d initial columns were added.\n", nColsAdded);
		DSPdebugMessage("master objective value %e, lower bound %e\n", primobj_, curlb);

		if (nColsAdded > 0) {
#if 1
			/** update master */
//			DSP_RTN_CHECK_RTN_CODE(updateModel(price, curlb));
			CoinZeroN(tr_center_, nrows_ - nrows_conv_);
			updateTrustRegion();
			dualobj_ = curlb;

//			/** use primal simplex after column generation */
//			if (!useCpxBarrier_)
//				si_->setHintParam(OsiDoDualInResolve, false);
//
//			/** resolve */
//			si_->resolve();
//
//			/** calculate primal objective value */
//			DSP_RTN_CHECK_RTN_CODE(calculatePrimalObjective());
//
//			double relgap = (primobj_-dualobj_)/(1.0e-10+fabs(primobj_))*100;
//			message_->print(2, "[Phase %d] Iteration %3d: Master objective %e, ", phase_, itercnt_, primobj_);
//			if (phase_ == 2)
//				message_->print(2, "Lb %e (gap %.2f %%), ", dualobj_, relgap);
#endif
			/** solve */
			DSP_RTN_CHECK_RTN_CODE(gutsOfSolve());
		} else {
			dualobj_ = primobj_;
		}
#if 0
		/** TODO: Why do we have this situation? Let's turn off the dual trust region for now. */
		if ((status_ == DSP_STAT_OPTIMAL || status_ == DSP_STAT_FEASIBLE)
				&& primobj_ - dualobj_ < -1.0e-6) {
			message_->print(1, "Numerical instability is detected. Turn off the dual TR device for this node.\n");
			/** turn off the TR device */
			for (int j = 0; j < ncols_tr_; ++j)
				si_->setColUpper(j, 0.0);

			/** reset dual objective */
			dualobj_ = -COIN_DBL_MAX;

			/** solve with TR */
			DSP_RTN_CHECK_RTN_CODE(gutsOfSolve());

			/** turn on the TR device */
			for (int j = 0; j < ncols_tr_; ++j)
				si_->setColUpper(j, COIN_DBL_MAX);
		}
#endif

	} else
		status_ = DSP_STAT_PRIM_INFEASIBLE;

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

DSP_RTN_CODE DwMasterTr::calculatePrimalObjective() {
	BGN_TRY_CATCH
#if 0
	if (phase_ == 1) {
		primobj_ = si_->getObjValue();
	} else {
		primobj_ = 0.0;
		for (int j = ncols_tr_; j < si_->getNumCols(); ++j)
			primobj_ += si_->getObjCoefficients()[j] * si_->getColSolution()[j];
	}
#else
	primobj_ = si_->getObjValue();
#endif
	/** get the sum of TR variable values. */
	tr_sumval_ = 0.0;
	for (int j = 0; j < ncols_tr_; ++j)
		tr_sumval_ += si_->getColSolution()[j];

#if 0
	for (int i = 0; i < nrows_orig_ + nrows_branch_; ++i)
		if (fabs(si_->getRowPrice()[i] - tr_center_[i]) >= tr_size_)
			printf("TR[%d]: |%e - %e| > [%e]\n", i, si_->getRowPrice()[i], tr_center_[i], tr_size_);
#endif
	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

DSP_RTN_CODE DwMasterTr::restoreCols() {
#define FREE_MEMORY \
	FREE_PTR(ws)

	CoinWarmStartBasis* ws = NULL;
	std::vector<int> existingCols;
	std::vector<int> delCols;

	BGN_TRY_CATCH

	/** resolve and get basis */
	si_->resolve();
	if (si_->isProvenOptimal())
		ws = dynamic_cast<CoinWarmStartBasis*>(si_->getWarmStart());

	/** mark the existing columns to match with basis; and mark columns to delete */
	for (unsigned k = 0, j = ncols_tr_; k < cols_generated_.size(); ++k)
		if (cols_generated_[k]->active_) {
			existingCols.push_back(k);
			delCols.push_back(j++);
		}

	if (delCols.size() > 0) {
		/** delete columns */
		si_->deleteCols(delCols.size(), &delCols[0]);

		/** add columns */
		for (unsigned k = 0; k < cols_generated_.size(); ++k) {
			cols_generated_[k]->active_ = true;
			si_->addCol(cols_generated_[k]->col_, cols_generated_[k]->lb_, cols_generated_[k]->ub_, cols_generated_[k]->obj_);
		}

		/** set warm start information */
		if (ws) {
			CoinWarmStartBasis* ews = dynamic_cast<CoinWarmStartBasis*>(si_->getEmptyWarmStart());
			for (int i = 0; i < ws->getNumArtificial(); ++i)
				ews->setArtifStatus(i, ws->getArtifStatus(i));
			for (int j = 0; j < ncols_tr_; ++j)
				ews->setStructStatus(j, ws->getStructStatus(j));
			for (unsigned k = 0; k < existingCols.size(); ++k)
				ews->setStructStatus(ncols_tr_ + existingCols[k], ws->getStructStatus(ncols_tr_+k));
			si_->setWarmStart(ews);
			FREE_PTR(ews)
		}
	}

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)
	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

DSP_RTN_CODE DwMasterTr::reduceCols() {
	BGN_TRY_CATCH

	if (phase_ == 2) {
		std::vector<int> delcols;
		for (unsigned k = 0, j = ncols_tr_; k < cols_generated_.size(); ++k) {
			if (cols_generated_[k]->active_) {
				/** age? */
				if (si_->getReducedCost()[j] < 1.0e-8)
					cols_generated_[k]->age_ = 0;
				else
					cols_generated_[k]->age_++;
				/** reduced cost fixing */
				if (cols_generated_[k]->age_ >= par_->getIntParam("DW/MASTER/COL_AGE_LIM") ||
						dualobj_ + si_->getReducedCost()[j] - bestprimobj_ > -1.0e-10) {
					cols_generated_[k]->active_ = false;
					delcols.push_back(j);
				}
				j++;
			}
		}
		if (delcols.size() > 0 && useBarrier_ == false) {
			CoinWarmStartBasis* ws = dynamic_cast<CoinWarmStartBasis*>(si_->getWarmStart());
			ws->deleteColumns(delcols.size(), &delcols[0]);
			si_->deleteCols(delcols.size(), &delcols[0]);
			si_->setWarmStart(ws);
			FREE_PTR(ws)
			//message_->print(1, "  Deleted %u columns.\n", delcols.size());

			si_->resolve();
			ws = dynamic_cast<CoinWarmStartBasis*>(si_->getWarmStart());
		}
	}
	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

bool DwMasterTr::terminationTest(int nnewcols, int itercnt, double relgap) {

	bool term = false;
	/** TODO: parameterize */
	double gaplim = 0.0001;

	/** time ticking */
	ticToc();

	if (!isTrBoundary()) {
		if (time_remains_ < 0 || iterlim_ <= itercnt) {
			status_ = DSP_STAT_LIM_ITERorTIME;
			term = true;
		} else if (nnewcols == 0 ||
				relgap < gaplim ||
				dualobj_ - bestprimobj_ > -1.0e-10) {
			status_ = DSP_STAT_FEASIBLE;
			term = true;
		}
	}

	if (!term && relgap < gaplim) {
		tr_size_ = CoinMin(2. * tr_size_, 1.0e+4);
		updateTrustRegion();
		message_->print(3, "Termination test: increased trust region size %e\n", tr_size_);
	}

	return term;
}

DSP_RTN_CODE DwMasterTr::updateModel(
		const double* price, /**< [in] price */
		double curLb         /**< [in] current lower bound */) {
	BGN_TRY_CATCH

	//DSPdebugMessage("current primobj %e\n", primobj_);

	if (curLb >= dualobj_ + 1.0e-4 * (primobj_ - dualobj_)) {
		message_->print(3, "  [TR] SERIOUS STEP: dual objective %e", curLb);

		/** increase trust region? */
		if (isTrBoundary() && curLb >= dualobj_ + 0.5 * (primobj_ - dualobj_)) {
			tr_size_ = CoinMin(2. * tr_size_, 1.0e+4);
			message_->print(3, ", increased trust region size %e", tr_size_);
		}

		/** update proximal point */
		CoinCopyN(price, nrows_orig_, tr_center_);
		CoinCopyN(price + nrows_core_, nrows_branch_, tr_center_ + nrows_orig_);
		message_->print(3, ", updated proximal point");

		/** update trust region */
		updateTrustRegion();

		/** update dual bound */
		dualobj_ = curLb;
		tr_cnt_ = 0;
		message_->print(3, "\n");
	} else {
		message_->print(3, "  [TR] null step: dual objective %e", curLb);
		if (isTrBoundary()) {
			/** increase trust region? */
			if (primobj_ < dualobj_) {
				tr_size_ = CoinMin(2.0 * tr_size_, 1.0e+4);
				message_->print(3, ", increased trust region size %e", tr_size_);
				/** update trust region */
				updateTrustRegion();
			}
		} else {
			double rho = CoinMin(1.0, tr_size_) * (dualobj_ - curLb) / (primobj_ - dualobj_);
			if (rho > 0) tr_cnt_++;
			if (rho >= 3 || (tr_cnt_ >= 3 && fabs(rho - 2.) < 1.0))
			{
				/** decrease trust region */
				tr_size_ *= 1.0 / CoinMin(rho, 4.);
				message_->print(3, ", decreased trust region size %e", tr_size_);
				tr_cnt_ = 0;

				/** update trust region */
				updateTrustRegion();
			}
		}
		message_->print(3, "\n");
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

DSP_RTN_CODE DwMasterTr::updateTrustRegion() {
	BGN_TRY_CATCH

	int half_ncols_tr = ncols_tr_ * 0.5;
	for (int j = 0; j < half_ncols_tr; ++j) {
		si_->setObjCoeff(j*2, tr_center_[j] + tr_size_);
		si_->setObjCoeff(j*2+1, -tr_center_[j] + tr_size_);
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}
