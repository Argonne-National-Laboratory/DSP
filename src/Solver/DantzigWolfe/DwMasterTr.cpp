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
#include "CoinUtility.hpp"
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
	rlbd_branch_ = new double [nrows_branch_];
	rubd_branch_ = new double [nrows_branch_];
	tr_center_ = new double [nrows_];
	starts_tr = new int [ncols_tr_ + 1];
	rows_tr = new int [ncols_tr_];
	elements_tr = new double [ncols_tr_];

	/** initial trust region center */
	CoinZeroN(tr_center_, nrows_orig_ + nrows_branch_);
	CoinFillN(tr_center_ + nrows_orig_ + nrows_branch_, nrows_conv_, COIN_DBL_MAX);

	/** create column-wise matrix and set number of rows */
	mat = new CoinPackedMatrix(true, 0, 0);
	mat->setDimensions(nrows_, 0);

	/** add variables related to trust region */
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
	for (int i = 0; i < nrows_branch_; ++i) {
		int j = branch_row_to_col_[nrows_orig_ + i];
		rlbd[nrows_orig_ + i] = org_clbd_[j];
		rubd[nrows_orig_ + i] = org_cubd_[j];
		rlbd_branch_[i] = org_clbd_[j];
		rubd_branch_[i] = org_cubd_[j];
	}
	CoinFillN(rlbd + nrows_orig_ + nrows_branch_, nrows_conv_, 1.0);
	CoinFillN(rubd + nrows_orig_ + nrows_branch_, nrows_conv_, 1.0);

	/** create solver */
	si_ = new OsiCpxSolverInterface();

	/** message setting */
	//dynamic_cast<OsiClpSolverInterface*>(si_)->getModelPtr()->setLogLevel(0);
	si_->messageHandler()->logLevel(0);
	//si_->passInMessageHandler(si_->messageHandler());
	//DSPdebug(si_->messageHandler()->logLevel(4));

	/** load problem data */
	si_->loadProblem(*mat, clbd, cubd, obj, rlbd, rubd);

	/** write mps */
	DSPdebug(si_->writeMps("master"));

	/** NOTE: This does not go to Phase 1. */
	phase_ = 2;

	/** set hint parameters */
	useCpxBarrier_ = false;

	cpx = dynamic_cast<OsiCpxSolverInterface*>(si_);
	if (cpx)
		CPXsetintparam(cpx->getEnvironmentPtr(), CPX_PARAM_THREADS, 1);

	if (useCpxBarrier_) {
		/** use barrier */
		CPXsetintparam(cpx->getEnvironmentPtr(), CPX_PARAM_LPMETHOD, CPX_ALG_BARRIER);
		/** no crossover */
		CPXsetintparam(cpx->getEnvironmentPtr(), CPX_PARAM_BARCROSSALG, -1);
		CPXsetintparam(cpx->getEnvironmentPtr(), CPX_PARAM_BARDISPLAY, 0);
	} else {
		si_->setHintParam(OsiDoPresolveInResolve, true);
		si_->setHintParam(OsiDoDualInResolve, false);
		si_->setHintParam(OsiDoInBranchAndCut, false);
	}

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

	DSP_RTN_CHECK_RTN_CODE(solvePhase1());

	if (primobj_ > 1.0e-8)
		status_ = DSP_STAT_PRIM_INFEASIBLE;
	else {

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
					/** do not consider inactive columns */
					if (cols_generated_[k]->active_ == false) continue;
					CoinPackedVector xlam = cols_generated_[k]->x_ * si_->getColSolution()[j];
					for (int i = 0; i < xlam.getNumElements(); ++i)
						primsol_[xlam.getIndices()[i]] += xlam.getElements()[i];
					j++;
				}

				/** heuristics */
				DSP_RTN_CHECK_RTN_CODE(heuristics());
				break;
			} else if (status_ == DSP_STAT_DUAL_INFEASIBLE)
				break;
		}
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DwMasterTr::solvePhase1() {
#define FREE_MEMORY \
	FREE_ARRAY_PTR(price) \
	FREE_ARRAY_PTR(piA)

	double* price = NULL;
	double* piA = NULL;

	BGN_TRY_CATCH

	phase_ = 1;
	status_ = DSP_STAT_FEASIBLE;

	/** set the objective function coefficient for phase 1 */
	for (int j = 0; j < ncols_tr_; ++j) {
		si_->setObjCoeff(j, 1.0);
		si_->setColUpper(j, COIN_DBL_MAX);
	}
	for (int j = ncols_tr_; j < si_->getNumCols(); ++j)
		si_->setObjCoeff(j, 0.0);

	if (ncols_tr_ == si_->getNumCols()) {
		DSPdebugMessage("Generate initial columns.\n");

		/** allocate memory */
		price = new double [nrows_];
		piA = new double [ncols_orig_];
		CoinZeroN(price, nrows_ - nrows_conv_);
		CoinFillN(price + nrows_ - nrows_conv_, nrows_conv_, COIN_DBL_MAX);

		/** generate columns */
		double dummy_lb;
		int nColsAdded = 0;
		DSP_RTN_CHECK_RTN_CODE(generateCols(price, piA, dummy_lb, nColsAdded));
		DSPdebugMessage("%d initial columns were added.\n", nColsAdded);
	}

	if (status_ == DSP_STAT_FEASIBLE) {
		/** solve */
		si_->resolve();
		if (si_->isProvenOptimal() == false)
			throw "Unexpected phase 1 status.";
		primobj_ = si_->getObjValue();
		DSPdebugMessage("Initial phase1 objective value is %e.\n", primobj_);

		if (primobj_ > 1.0e-10) {
			/** do column generation for phase 1 */
			DSP_RTN_CHECK_RTN_CODE(DwAlgo::gutsOfSolve());
		}
	} else
		primobj_ = COIN_DBL_MAX;

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

DSP_RTN_CODE DwMasterTr::solvePhase2() {
#define FREE_MEMORY \
	FREE_ARRAY_PTR(price) \
	FREE_ARRAY_PTR(piA) \
	FREE_PTR(ws)

	double* price = NULL;
	double* piA = NULL;
	int nColsAdded = 0;

	std::vector<double> prevsol;
	CoinWarmStartBasis* ws = NULL;

	BGN_TRY_CATCH

	phase_ = 2;

	/** set the objective function coefficient for phase 2 */
	for (int j = 0; j < ncols_tr_; ++j) {
		si_->setObjCoeff(j, 0.0);
		si_->setColUpper(j, 0.0);
	}
	for (unsigned k = 0, j = ncols_tr_; k < cols_generated_.size(); ++k)
		if (cols_generated_[k]->active_)
			si_->setObjCoeff(j++, cols_generated_[k]->obj_);

	/** allocate memory */
	price = new double [nrows_];
	piA = new double [ncols_orig_];
	prevsol.reserve(si_->getNumCols());

	/** solve without the dual trust-region */
	si_->resolve();
	if (si_->isProvenOptimal() == false)
		throw "Unexpected phase 2 status.";

	/** calculate primal objective value */
	DSP_RTN_CHECK_RTN_CODE(calculatePrimalObjective());

	CoinCopyN(si_->getRowPrice(), nrows_, price);

	/** get warm-start information */
	if (!useCpxBarrier_)
		DSP_RTN_CHECK_RTN_CODE(getWarmStartInfo(prevsol, ws));

	/** generate columns */
	DSP_RTN_CHECK_RTN_CODE(generateCols(price, piA, dualobj_, nColsAdded));

	/** termination test: subproblem solution may decalre infeasible. */
	if (status_ == DSP_STAT_OPTIMAL || status_ == DSP_STAT_FEASIBLE) {
		DSPdebugMessage("%d initial columns were added.\n", nColsAdded);
		DSPdebugMessage("master objective value %e, lower bound %e\n", primobj_, dualobj_);

		if (nColsAdded > 0) {
			/** put the trust-region variables back */
			for (int j = 0; j < ncols_tr_; ++j)
				si_->setColUpper(j, COIN_DBL_MAX);

			/** set trust region */
			tr_size_ = 10.0;
			CoinCopyN(price, nrows_, tr_center_);
			DSP_RTN_CHECK_RTN_CODE(updateTrustRegion());

			/** set warm-start information */
			if (!useCpxBarrier_)
				DSP_RTN_CHECK_RTN_CODE(setWarmStartInfo(prevsol, ws));

			/** solve */
			DSP_RTN_CHECK_RTN_CODE(DwAlgo::gutsOfSolve());
		} else {
			dualobj_ = primobj_;
		}

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
			DSP_RTN_CHECK_RTN_CODE(DwAlgo::gutsOfSolve());

			/** turn on the TR device */
			for (int j = 0; j < ncols_tr_; ++j)
				si_->setColUpper(j, COIN_DBL_MAX);
		}
	} else
		status_ = DSP_STAT_PRIM_INFEASIBLE;

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

DSP_RTN_CODE DwMasterTr::calculatePrimalObjective() {
	BGN_TRY_CATCH
#if 1
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
				/** reduced cost fixing */
				if (dualobj_ + si_->getReducedCost()[j] - bestprimobj_ > -1.0e-10) {
					cols_generated_[k]->active_ = false;
					delcols.push_back(j);
				}
				j++;
			}
		}
		if (delcols.size() > 0) {
			CoinWarmStartBasis* ws = dynamic_cast<CoinWarmStartBasis*>(si_->getWarmStart());
			ws->deleteColumns(delcols.size(), &delcols[0]);
			si_->deleteCols(delcols.size(), &delcols[0]);
			si_->setWarmStart(ws);
			FREE_PTR(ws)
			//message_->print(1, "  Deleted %u columns.\n", delcols.size());
		}
	}
	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

bool DwMasterTr::terminationTest(int nnewcols, int itercnt, double relgap) {

	bool term = false;
	int heuristicFreq = 10;

	/** time ticking */
	ticToc();

	if (phase_ == 1) {
		if (time_remains_ < 0 || iterlim_ <= itercnt) {
			status_ = DSP_STAT_LIM_ITERorTIME;
			term = true;
		} else if (si_->isProvenOptimal() && primobj_ < 1.0e-8) {
			status_ = DSP_STAT_FEASIBLE;
			DSPdebugMessage("Phase 1 found a feasible solution! %e < %e\n", si_->getObjValue(), 1.0e-8);
			term = true;
		}
	} else {
		if (time_remains_ < 0 || iterlim_ <= itercnt) {
			status_ = DSP_STAT_LIM_ITERorTIME;
			term = true;
		} else if (relgap < 0.01 - 1.0e-8 ||
			(phase_ == 2 && dualobj_ - bestprimobj_ > -1.0e-10)) {
			status_ = DSP_STAT_FEASIBLE;
			term = true;
		}
	}

	return term;
}

DSP_RTN_CODE DwMasterTr::updateModel(
		const double* price, /**< [in] price */
		double curLb         /**< [in] current lower bound */) {
	BGN_TRY_CATCH

	//DSPdebugMessage("current primobj %e\n", primobj_);

	if (phase_ == 1)
		DwAlgo::updateModel(price, curLb);
	else {
		if (curLb >= dualobj_ + 1.0e-4 * (primobj_ - dualobj_)) {
			message_->print(3, "  [TR] SERIOUS STEP: dual objective %e", curLb);

			if (dualobj_ > -COIN_DBL_MAX) {
				/** increase trust region? */
				if (isTrBoundary(price) && curLb >= dualobj_ + 0.5 * (primobj_ - dualobj_)) {
					tr_size_ = CoinMin(2. * tr_size_, 1.0e+4);
					message_->print(3, ", increased trust region size %e", tr_size_);
				}

				/** update proximal point */
				CoinCopyN(price, nrows_, tr_center_);
				message_->print(3, ", updated proximal point");

				/** update trust region */
				updateTrustRegion();
			}

			/** update dual bound */
			dualobj_ = curLb;
			tr_cnt_ = 0;
			message_->print(3, "\n");
		} else {
			message_->print(3, "  [TR] null step: dual objective %e", curLb);
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
			message_->print(3, "\n");
		}
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

bool DwMasterTr::isTrBoundary(const double* price) {
	bool isBoundary = false;

	BGN_TRY_CATCH

	for (int j = 0; j < nrows_orig_ + nrows_branch_; ++j) {
		if (price[2*j] + price[2*j+1] > 1.0e-10) {
			isBoundary = true;
			break;
		}
	}

	END_TRY_CATCH_RTN(;,isBoundary)

	return isBoundary;
}

DSP_RTN_CODE DwMasterTr::updateTrustRegion() {
	BGN_TRY_CATCH

	for (int j = 0; j < nrows_orig_ + nrows_branch_; ++j) {
		si_->setObjCoeff(j*2, tr_center_[j] + tr_size_);
		si_->setObjCoeff(j*2+1, -tr_center_[j] + tr_size_);
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}
