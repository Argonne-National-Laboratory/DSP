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

	/** allocate memory */
	obj = new double [ncols_tr_];
	clbd = new double [ncols_tr_];
	cubd = new double [ncols_tr_];
	rlbd = new double [nrows_branch_];
	rubd = new double [nrows_branch_];
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

DSP_RTN_CODE DwMasterTr::heuristics() {
	if (!runHeuristics_)
		return DSP_RTN_OK;

	double stime;

	BGN_TRY_CATCH

	if (iterlim_ == 1 && par_->getBoolParam("DW/HEURISTICS/TRIVIAL")) {
		message_->print(1, "Heuristic (trivial) searches solutions...\n");
		stime = CoinGetTimeOfDay();
		DSP_RTN_CHECK_RTN_CODE(heuristicTrivial());
		message_->print(1, "Heuristic (trivial) spent %.2f seconds [best %e].\n", CoinGetTimeOfDay() - stime, bestprimobj_);
	}

//	DSPdebugMessage("Running FP-like heuristic with direction 1.\n");
//	stime = CoinGetTimeOfDay();
//	DSP_RTN_CHECK_RTN_CODE(heuristicFp(1));
//	message_->print(1, "Heuristic (FP-like[+1]) spent %.2f seconds [best %e].\n", CoinGetTimeOfDay() - stime, bestprimobj_);
//
//	DSPdebugMessage("Running FP-like heuristic with direction -1.\n");
//	stime = CoinGetTimeOfDay();
//	DSP_RTN_CHECK_RTN_CODE(heuristicFp(-1));
//	message_->print(1, "Heuristic (FP-like[-1]) spent %.2f seconds [best %e].\n", CoinGetTimeOfDay() - stime, bestprimobj_);

	if (par_->getBoolParam("DW/HEURISTICS/DIVE")) {
		message_->print(1, "Heuristic (dive) searches solutions...\n");
		stime = CoinGetTimeOfDay();
		DSP_RTN_CHECK_RTN_CODE(heuristicDive());
		message_->print(1, "Heuristic (dive) spent %.2f seconds [best %e].\n", CoinGetTimeOfDay() - stime, bestprimobj_);
	}

	/** restore the original settings */
	iterlim_ = par_->getIntParam("DW/ITER_LIM");

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

DSP_RTN_CODE DwMasterTr::heuristicFp(int direction) {
#define FREE_MEMORY \
	FREE_ARRAY_PTR(rlbd) \
	FREE_ARRAY_PTR(rubd) \
	FREE_ARRAY_PTR(primsol_org)

	/** row bounds for branching rows */
	double* rlbd = NULL;
	double* rubd = NULL;

	/** original solutions */
	double* primsol_org = NULL;
	double primobj_org, dualobj_org;
	int status_org;

	BGN_TRY_CATCH

	DSP_RTN_CHECK_RTN_CODE(
			preHeuristic(rlbd, rubd, primsol_org, primobj_org, dualobj_org, status_org));

	/** round and fix bounds */
	for (int i = 0; i < nrows_branch_; ++i) {
		/** do not consider those with non-positive objective coefficients */
		int colind = branch_row_to_col_[nrows_orig_ + i];
		if (org_obj_[colind] * direction <= 0) continue;

		double rounded = round(primsol_[colind]);
		si_->setRowBounds(nrows_orig_ + i, rounded, rounded);
		worker_->setColBounds(1, &colind, &rounded, &rounded);
	}

	/** solve */
	DSP_RTN_CHECK_RTN_CODE(solvePhase1());

	if (primobj_ > 1.0e-8)
		status_ = DSP_STAT_PRIM_INFEASIBLE;
	else
		DSP_RTN_CHECK_RTN_CODE(solvePhase2());

	/** collect solutions */
	if (status_ == DSP_STAT_OPTIMAL || status_ == DSP_STAT_FEASIBLE) {
		/** check integer feasibility */
		bool fractional = false;
		for (int i = 0; i < nrows_branch_; ++i) {
			double x = si_->getRowActivity()[nrows_orig_ + i];
			if (fabs(x - floor(x + 0.5)) > 1.0e-10) {
				fractional = true;
				DSPdebugMessage("Heuristic found a fractional solution.\n");
				break;
			}
		}
		if (!fractional) {
			DSPdebugMessage("Heuristic found an integer solution.\n");
			if (bestprimobj_ > primobj_) {
				bestprimobj_ = primobj_;
				/** recover original solution */
				CoinZeroN(bestprimsol_, ncols_orig_);
				for (unsigned k = 0, j = ncols_tr_; k < cols_generated_.size(); ++k) {
					/** do not consider inactive columns */
					if (cols_generated_[k]->active_ == false)
						continue;
					CoinPackedVector xlam = cols_generated_[k]->x_ * si_->getColSolution()[j];
					for (int i = 0; i < xlam.getNumElements(); ++i)
						bestprimsol_[xlam.getIndices()[i]] += xlam.getElements()[i];
					j++;
				}
				DSPdebugMessage("Heuristic updated the best upper bound %e.\n", bestprimobj_);
			}
		}
	} else if (status_ == DSP_STAT_DUAL_INFEASIBLE) {
		//DSPdebug(si_->writeMps("master"));
	}

	DSP_RTN_CHECK_RTN_CODE(
			postHeuristic(rlbd, rubd, primsol_org, primobj_org, dualobj_org, status_org));

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

DSP_RTN_CODE DwMasterTr::heuristicTrivial() {
#define FREE_MEMORY \
	FREE_ARRAY_PTR(rlbd) \
	FREE_ARRAY_PTR(rubd) \
	FREE_ARRAY_PTR(primsol_org)

	/** row bounds for branching rows */
	double* rlbd = NULL;
	double* rubd = NULL;

	/** original solutions */
	double* primsol_org = NULL;
	double primobj_org, dualobj_org;
	int status_org;
	int iterlim;

	BGN_TRY_CATCH

	DSP_RTN_CHECK_RTN_CODE(
			preHeuristic(rlbd, rubd, primsol_org, primobj_org, dualobj_org, status_org));

	/** fix bounds */
	std::vector<int> branchIndices;
	std::vector<double> branchBounds;
	for (int i = 0; i < nrows_branch_; ++i) {
		/** do not consider those with negative objective coefficients */
		if (org_obj_[branch_row_to_col_[nrows_orig_ + i]] >= 0) continue;
		/** skip fixed bounds */
		if (rlbd[i] == rubd[i]) continue;
		si_->setRowBounds(nrows_orig_ + i, 0.0, 0.0);
		branchIndices.push_back(branch_row_to_col_[nrows_orig_ + i]);
		branchBounds.push_back(0.0);
	}
	worker_->setColBounds(branchIndices.size(), &branchIndices[0], &branchBounds[0], &branchBounds[0]);

	/** solve */
	itercnt_ = 0;
	iterlim_ = par_->getIntParam("DW/HEURISTICS/TRIVIAL/ITER_LIM");
	time_remains_ = par_->getDblParam("DW/HEURISTICS/TRIVIAL/TIME_LIM");
	DSP_RTN_CHECK_RTN_CODE(solvePhase1());

	if (primobj_ > 1.0e-8)
		status_ = DSP_STAT_PRIM_INFEASIBLE;
	else {
		itercnt_ = 0;
		DSP_RTN_CHECK_RTN_CODE(solvePhase2());
	}

	/** collect solutions */
	if (status_ == DSP_STAT_OPTIMAL ||
			status_ == DSP_STAT_FEASIBLE ||
			status_ == DSP_STAT_LIM_ITERorTIME) {
		/** check integer feasibility */
		bool fractional = false;
		for (int i = 0; i < nrows_branch_; ++i) {
			double x = si_->getRowActivity()[nrows_orig_ + i];
			if (fabs(x - floor(x + 0.5)) > 1.0e-10) {
				fractional = true;
				DSPdebugMessage("Heuristic found a fractional solution.\n");
				break;
			}
		}
		if (!fractional) {
			message_->print(1, "Heuristic found an integer solution (objective %e).\n", primobj_);
			if (bestprimobj_ > primobj_) {
				bestprimobj_ = primobj_;
				/** recover original solution */
				CoinZeroN(bestprimsol_, ncols_orig_);
				for (unsigned k = 0, j = ncols_tr_; k < cols_generated_.size(); ++k) {
					/** do not consider inactive columns */
					if (cols_generated_[k]->active_ == false)
						continue;
					CoinPackedVector xlam = cols_generated_[k]->x_ * si_->getColSolution()[j];
					for (int i = 0; i < xlam.getNumElements(); ++i)
						bestprimsol_[xlam.getIndices()[i]] += xlam.getElements()[i];
					j++;
				}
				message_->print(1, "Heuristic updated the best upper bound %e.\n", bestprimobj_);
			}
		}
	} else if (status_ == DSP_STAT_DUAL_INFEASIBLE) {
		//DSPdebug(si_->writeMps("master"));
	}

	DSP_RTN_CHECK_RTN_CODE(
			postHeuristic(rlbd, rubd, primsol_org, primobj_org, dualobj_org, status_org));

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

DSP_RTN_CODE DwMasterTr::heuristicDive() {
#define FREE_MEMORY \
	FREE_ARRAY_PTR(rlbd) \
	FREE_ARRAY_PTR(rubd) \
	FREE_ARRAY_PTR(primsol_org)

	/** row bounds for branching rows */
	double* rlbd = NULL;
	double* rubd = NULL;

	/** original solutions */
	double* primsol_org = NULL;
	double primobj_org, dualobj_org;
	int status_org;

	std::vector<CoinTriple<int,int,double> > branchList;

	BGN_TRY_CATCH

	DSP_RTN_CHECK_RTN_CODE(
			preHeuristic(rlbd, rubd, primsol_org, primobj_org, dualobj_org, status_org));

	DSP_RTN_CHECK_RTN_CODE(gutsOfDive(branchList, 0));

	DSP_RTN_CHECK_RTN_CODE(
			postHeuristic(rlbd, rubd, primsol_org, primobj_org, dualobj_org, status_org));

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

DSP_RTN_CODE DwMasterTr::gutsOfDive(
		std::vector<CoinTriple<int,int,double> > branchList,
		int depth) {

	/** parameters for dive with backtracking */
	int maxdiscrepancy = 2;
	int maxdepth = 3;

	int findPhase = 0;
	int branchIndex, branchDirection, solvePhase;
	double branchValue;
	double dist, mindist;
	bool isInteger;

	int status; /**< status at the current depth */

	BGN_TRY_CATCH

	while (1) {

		if (depth > 0) {
			/** solve */
			itercnt_ = 0;
			iterlim_ = par_->getIntParam("DW/HEURISTICS/DIVE/ITER_LIM");
			time_remains_ = par_->getDblParam("DW/HEURISTICS/DIVE/TIME_LIM");
			DSP_RTN_CHECK_RTN_CODE(solvePhase1());

			/** recover original solution */
			CoinZeroN(primsol_, ncols_orig_);
			for (unsigned k = 0, j = ncols_tr_; k < cols_generated_.size(); ++k) {
				/** do not consider inactive columns */
				if (cols_generated_[k]->active_ == false)
					continue;
				CoinPackedVector xlam = cols_generated_[k]->x_ * si_->getColSolution()[j];
				for (int i = 0; i < xlam.getNumElements(); ++i)
					primsol_[xlam.getIndices()[i]] += xlam.getElements()[i];
				j++;
			}
		}

		/** find a fractional */
		findPhase = 0;
		branchIndex = -1;
		mindist = 1.0;
		isInteger = true;
		while (findPhase < 2 && branchIndex < 0) {
			findPhase++;
			for (int i = 0; i < nrows_branch_; ++i) {
				int colind = branch_row_to_col_[nrows_orig_ + i];
				/** do not consider those with non-negative objective coefficients */
				if (findPhase == 0 && org_obj_[colind] >= 0) continue;
				/** do not consider those with negative objective coefficients */
				if (findPhase == 1 && org_obj_[colind] < 0) continue;

				/** skip if already fixed */
				if (si_->getRowUpper()[nrows_orig_ + i] - si_->getRowLower()[nrows_orig_ + i] < 1.0e-8)
					continue;

				dist = fabs(primsol_[colind] - floor(primsol_[colind] + 0.5));
				if (dist < 1.0e-8) continue;

				/** mark as not integer */
				isInteger = false;

				double candValue = primsol_[colind];
				double candRounded = round(candValue);
				int candDirection = (candRounded < candValue) ? -1 : 1;

				/** check if the branch found is in branchList. */
				for (unsigned j = 0; j < branchList.size(); ++j) {
					if (branchList[j].first == i &&
						branchList[j].second == candDirection &&
						fabs(branchList[j].third - candRounded) < 1.0e-10) {
#if 0
						/** yes, it exists; so do not branch on this. */
						dist = 1.0;
						break;
#else
						/** flip */
						if (candValue > candRounded)
							candRounded += 1.0;
						else
							candRounded -= 1.0;
						dist = 0.5 - dist;
						if (fabs(branchList[j].third - candRounded) < 1.0e-10) {
							dist = 1.0;
							break;
						}
#endif
					}
				}

				if (dist < mindist) {
					mindist = dist;
					branchDirection = candDirection;
					branchIndex = i;
					branchValue = candValue;
				}
			}
		}

		/** found a fractional variable */
		if (branchIndex > -1) {
			/** keep the current node bounds */
			double rlbd_node = si_->getRowLower()[nrows_orig_ + branchIndex];
			double rubd_node = si_->getRowUpper()[nrows_orig_ + branchIndex];

			/** fix bounds */
			double rounded = round(branchValue);
			si_->setRowBounds(nrows_orig_ + branchIndex, rounded, rounded);
			worker_->setColBounds(1, &branch_row_to_col_[nrows_orig_ + branchIndex], &rounded, &rounded);
			message_->print(2, "Diving fixed variable %d [%e] to %e (discrepancy %u, depth %d).\n", branchIndex, branchValue, rounded, branchList.size(), depth);

			/** recursive call */
			status = status_;
			DSP_RTN_CHECK_RTN_CODE(gutsOfDive(branchList, depth+1));
			status_ = status;

			/** restore node bounds */
			si_->setRowBounds(nrows_orig_ + branchIndex, rlbd_node, rubd_node);
			worker_->setColBounds(1, &branch_row_to_col_[nrows_orig_ + branchIndex], &rlbd_node, &rubd_node);

			/** put a branch to the list */
			branchList.push_back(CoinMakeTriple(branchIndex, branchDirection, branchValue));

			DSPdebugMessage("discrepancy %u, depth %d\n", branchList.size(), depth);
			if (branchList.size() > maxdiscrepancy || depth > maxdepth)
				break;
		} else if (!isInteger) {
			break;
		} else {
			message_->print(1, "Diving found an integer solution.\n");
#define FIX_NOW
#ifdef FIX_NOW
			/** backup and fix all bounds */
			std::vector<int> branchIndices;
			std::vector<double> branchBounds;
			double* rlbd_tmp = new double [nrows_branch_];
			double* rubd_tmp = new double [nrows_branch_];
			CoinCopyN(si_->getRowLower() + nrows_orig_, nrows_branch_, rlbd_tmp);
			CoinCopyN(si_->getRowUpper() + nrows_orig_, nrows_branch_, rubd_tmp);
			for (int j = 0; j < nrows_branch_; ++j) {
				double rounded = round(primsol_[branch_row_to_col_[nrows_orig_ + j]]);
				si_->setRowBounds(nrows_orig_ + j, rounded, rounded);
				branchIndices.push_back(branch_row_to_col_[nrows_orig_ + j]);
				branchBounds.push_back(rounded);
			}
			worker_->setColBounds(branchIndices.size(), &branchIndices[0], &branchBounds[0], &branchBounds[0]);

			/** solve */
			itercnt_ = 0;
			iterlim_ = par_->getIntParam("DW/HEURISTICS/DIVE/ITER_LIM");
			time_remains_ = par_->getDblParam("DW/HEURISTICS/DIVE/TIME_LIM");
			DSP_RTN_CHECK_RTN_CODE(solvePhase1());

			/** determine solution status */
			if (primobj_ > 1.0e-8) {
				message_->print(1, "The integer solution is infeasible.\n");
				status_ = DSP_STAT_PRIM_INFEASIBLE;
				break;
			} else
				status_ = DSP_STAT_FEASIBLE;
#endif
			message_->print(1, "Diving is evaluating the integer solution.\n");
			itercnt_ = 0;
			iterlim_ = par_->getIntParam("DW/HEURISTICS/DIVE/ITER_LIM");
			time_remains_ = par_->getDblParam("DW/HEURISTICS/DIVE/TIME_LIM");
			DSP_RTN_CHECK_RTN_CODE(solvePhase2());

			/** collect solutions */
			bool terminateLoop = false;
			if (status_ == DSP_STAT_OPTIMAL || status_ == DSP_STAT_FEASIBLE) {
				/** recover original solution */
				CoinZeroN(primsol_, ncols_orig_);
				for (unsigned k = 0, j = ncols_tr_; k < cols_generated_.size(); ++k) {
					/** do not consider inactive columns */
					if (cols_generated_[k]->active_ == false)
						continue;
					CoinPackedVector xlam = cols_generated_[k]->x_ * si_->getColSolution()[j];
					for (int i = 0; i < xlam.getNumElements(); ++i)
						primsol_[xlam.getIndices()[i]] += xlam.getElements()[i];
					j++;
				}
				bool fractional = false;
#ifndef FIX_NOW
				/** check integer feasibility */
				for (int i = 0; i < nrows_branch_; ++i) {
					double x = primsol_[branch_row_to_col_[nrows_orig_ + i]];
					if (fabs(x - floor(x + 0.5)) > 1.0e-10) {
						fractional = true;
						DSPdebugMessage("Heuristic found a fractional solution (x %d [%e]).\n", i, x);
						break;
					}
				}
#endif
				if (!fractional) {
					primobj_ = 0.0;
					for (int j = ncols_tr_; j < si_->getNumCols(); ++j)
						primobj_ += si_->getObjCoefficients()[j] * si_->getColSolution()[j];
					message_->print(1, "Diving found an integer solution (objective %e).\n", primobj_);
					if (bestprimobj_ > primobj_) {
						bestprimobj_ = primobj_;
						CoinCopyN(primsol_, ncols_orig_, bestprimsol_);
						message_->print(1, "Diving updated the best upper bound %e.\n", bestprimobj_);
					}
					terminateLoop = true;
				}
			}

#ifdef FIX_NOW
			for (int j = 0; j < nrows_branch_; ++j)
				si_->setRowBounds(nrows_orig_ + j, rlbd_tmp[j], rubd_tmp[j]);
			worker_->setColBounds(branchIndices.size(), &branchIndices[0], rlbd_tmp, rubd_tmp);
			FREE_ARRAY_PTR(rlbd_tmp)
			FREE_ARRAY_PTR(rubd_tmp)
#endif

			if (terminateLoop)
				break;
		}
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DwMasterTr::preHeuristic(
		double*& rlbd,    /**< [out] original row lower bounds */
		double*& rubd,    /**< [out] original row lower bounds */
		double*& primsol, /**< [out] original primal solution */
		double& primobj,  /**< [out] original primal objective */
		double& dualobj,  /**< [out] original dual objective */
		int& status       /**< [out] original solution status */) {
	BGN_TRY_CATCH

#ifndef DSP_DEBUG
	message_->logLevel_ = 1;
#endif

	/** save the original row bounds */
	rlbd = new double [nrows_branch_];
	rubd = new double [nrows_branch_];
	CoinCopyN(si_->getRowLower() + nrows_orig_, nrows_branch_, rlbd);
	CoinCopyN(si_->getRowUpper() + nrows_orig_, nrows_branch_, rubd);

	/** save the original solutions */
	primsol = new double [ncols_orig_];
	CoinCopyN(primsol_, ncols_orig_, primsol);
	primobj = primobj_;
	dualobj = dualobj_;
	status = status_;

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DwMasterTr::postHeuristic(
		double*& rlbd,    /**< [out] original row lower bounds */
		double*& rubd,    /**< [out] original row lower bounds */
		double*& primsol, /**< [out] original primal solution */
		double& primobj,  /**< [out] original primal objective */
		double& dualobj,  /**< [out] original dual objective */
		int& status       /**< [out] original solution status */) {
	BGN_TRY_CATCH

	message_->logLevel_ = par_->getIntParam("LOG_LEVEL");

	/** restore the original solutions */
	CoinCopyN(primsol, ncols_orig_, primsol_);
	primobj_ = primobj;
	dualobj_ = dualobj;
	status_ = status;

	/** restore the original row bounds */
	std::vector<int> branchIndices;
	for (int i = 0; i < nrows_branch_; ++i) {
		si_->setRowBounds(nrows_orig_ + i, rlbd[i], rubd[i]);
		branchIndices.push_back(branch_row_to_col_[nrows_orig_ + i]);
	}
	worker_->setColBounds(nrows_branch_, &branchIndices[0], rlbd, rubd);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
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
