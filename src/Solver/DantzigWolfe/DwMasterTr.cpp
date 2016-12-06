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
#define FREE_MEMORY       \
	FREE_PTR(mat)

	OsiCpxSolverInterface* cpx = NULL;

	/** master problem */
	CoinPackedMatrix * mat = NULL;
	std::vector<double> obj;
	std::vector<double> clbd;
	std::vector<double> cubd;
	std::vector<double> rlbd;
	std::vector<double> rubd;

	BGN_TRY_CATCH

	/** number of initial columns */
	ncols_tr_ = 2 * (nrows_orig_ + nrows_branch_);
	int ncols = cols_generated_.size();

	/** allocate memory */
	obj.reserve(ncols_tr_ + ncols);
	clbd.reserve(ncols_tr_ + ncols);
	cubd.reserve(ncols_tr_ + ncols);
	rlbd.reserve(nrows_branch_);
	rubd.reserve(nrows_branch_);
	rlbd_branch_ = new double [nrows_branch_];
	rubd_branch_ = new double [nrows_branch_];
	tr_center_ = new double [nrows_];

	/** initial trust region center */
	CoinZeroN(tr_center_, nrows_orig_ + nrows_branch_);
	CoinFillN(tr_center_ + nrows_orig_ + nrows_branch_, nrows_conv_, COIN_DBL_MAX);

	/** create column-wise matrix and set number of rows */
	mat = new CoinPackedMatrix(true, 0, 0);
	mat->setDimensions(nrows_, 0);

	/** add variables related to trust region */
	for (int j = 0; j < nrows_orig_ + nrows_branch_; ++j) {
		clbd.push_back(0.0);
		cubd.push_back(COIN_DBL_MAX);
		obj.push_back(tr_size_);
		double vecelem = 1.0;
		mat->appendCol(1, &j, &vecelem);

		clbd.push_back(0.0);
		cubd.push_back(COIN_DBL_MAX);
		obj.push_back(tr_size_);
		vecelem = -1.0;
		mat->appendCol(1, &j, &vecelem);
	}

	/** add initial columns */
	for (int j = 0; j < ncols; ++j) {
		clbd.push_back(cols_generated_[j]->lb_);
		cubd.push_back(cols_generated_[j]->ub_);
		obj.push_back(cols_generated_[j]->obj_);
		mat->appendCol(cols_generated_[j]->col_);
	}
	DSPdebug(mat->verifyMtx(4));

	/** Set row bounds */
	for (int i = 0; i < nrows_orig_; ++i) {
		rlbd.push_back(org_rlbd_[i]);
		rubd.push_back(org_rubd_[i]);
	}
	for (int i = 0; i < nrows_branch_; ++i) {
		int j = branch_row_to_col_[nrows_orig_ + i];
		rlbd.push_back(org_clbd_[j]);
		rubd.push_back(org_cubd_[j]);
		rlbd_branch_[i] = org_clbd_[j];
		rubd_branch_[i] = org_cubd_[j];
	}
	for (int i = 0; i < nrows_conv_; ++i) {
		rlbd.push_back(1.0);
		rubd.push_back(1.0);
	}

	/** create solver */
	si_ = new OsiCpxSolverInterface();

	/** message setting */
	//dynamic_cast<OsiClpSolverInterface*>(si_)->getModelPtr()->setLogLevel(0);
	si_->messageHandler()->logLevel(0);
	//si_->passInMessageHandler(si_->messageHandler());
	//DSPdebug(si_->messageHandler()->logLevel(4));

	/** load problem data */
	si_->loadProblem(*mat, &clbd[0], &cubd[0], &obj[0], &rlbd[0], &rubd[0]);

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
	BGN_TRY_CATCH

	DSP_RTN_CHECK_RTN_CODE(solvePhase1());

	if (primobj_ > 1.0e-8)
		status_ = DSP_STAT_PRIM_INFEASIBLE;
	else
		DSP_RTN_CHECK_RTN_CODE(solvePhase2());

	/** collect solutions */
	if (status_ == DSP_STAT_OPTIMAL || status_ == DSP_STAT_FEASIBLE) {
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
		if (runHeuristics_) {
			DSPdebugMessage("Running a trivial heuristic.\n");
			double stime = CoinGetTimeOfDay();
			DSP_RTN_CHECK_RTN_CODE(heuristicTrivial());
			message_->print(1, "Heuristic (trivial) spent %.2f seconds [best %e].\n", CoinGetTimeOfDay() - stime, bestprimobj_);

//			DSPdebugMessage("Running FP-like heuristic with direction 1.\n");
//			stime = CoinGetTimeOfDay();
//			DSP_RTN_CHECK_RTN_CODE(heuristicFp(1));
//			message_->print(1, "Heuristic (FP-like[+1]) spent %.2f seconds [best %e].\n", CoinGetTimeOfDay() - stime, bestprimobj_);
//
//			DSPdebugMessage("Running FP-like heuristic with direction -1.\n");
//			stime = CoinGetTimeOfDay();
//			DSP_RTN_CHECK_RTN_CODE(heuristicFp(-1));
//			message_->print(1, "Heuristic (FP-like[-1]) spent %.2f seconds [best %e].\n", CoinGetTimeOfDay() - stime, bestprimobj_);

			stime = CoinGetTimeOfDay();
			DSP_RTN_CHECK_RTN_CODE(heuristicDive());
			message_->print(1, "Heuristic (dive) spent %.2f seconds [best %e].\n", CoinGetTimeOfDay() - stime, bestprimobj_);
		}
	} else if (status_ == DSP_STAT_DUAL_INFEASIBLE) {
		//DSPdebug(si_->writeMps("master"));
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

	/** number of auxiliary columns for trust region */
	int ncols_tr = 2*(nrows_orig_+nrows_branch_);

	/** set the objective function coefficient for phase 1 */
	for (int j = 0; j < ncols_tr; ++j) {
		si_->setObjCoeff(j, 1.0);
		si_->setColUpper(j, COIN_DBL_MAX);
	}
	for (int j = ncols_tr; j < si_->getNumCols(); ++j)
		si_->setObjCoeff(j, 0.0);

	if (ncols_tr == si_->getNumCols()) {
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
	FREE_ARRAY_PTR(piA)

	double* price = NULL;
	double* piA = NULL;
	int nColsAdded = 0;

	BGN_TRY_CATCH

	phase_ = 2;

	/** number of auxiliary columns for trust region */
	int ncols_tr = 2*(nrows_orig_+nrows_branch_);

	/** allocate memory */
	price = new double [nrows_];
	piA = new double [ncols_orig_];

	/** set the objective function coefficient for phase 2 */
	for (int j = 0; j < ncols_tr; ++j) {
		si_->setObjCoeff(j, 0.0);
		si_->setColUpper(j, 0.0);
	}
	for (unsigned k = 0, j = ncols_tr; k < cols_generated_.size(); ++k)
		if (cols_generated_[k]->active_) {
			si_->setObjCoeff(j++, cols_generated_[k]->obj_);
		}

	/** solve without the dual trust-region */
	si_->resolve();
	if (si_->isProvenOptimal() == false)
		throw "Unexpected phase 2 status.";
	primobj_ = si_->getObjValue();
	CoinCopyN(si_->getRowPrice(), nrows_, price);

	/** generate columns */
	DSP_RTN_CHECK_RTN_CODE(generateCols(price, piA, dualobj_, nColsAdded));

	/** termination test: subproblem solution may decalre infeasible. */
	if (status_ == DSP_STAT_OPTIMAL || status_ == DSP_STAT_FEASIBLE) {
		DSPdebugMessage("%d initial columns were added.\n", nColsAdded);
		DSPdebugMessage("master objective value %e, lower bound %e\n", primobj_, dualobj_);
		if (nColsAdded > 0) {
			/** put the trust-region variables back */
			for (int j = 0; j < ncols_tr; ++j)
				si_->setColUpper(j, COIN_DBL_MAX);

			/** set trust region */
			tr_size_ = 10.0;
			CoinCopyN(price, nrows_, tr_center_);
			DSP_RTN_CHECK_RTN_CODE(updateTrustRegion());

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

bool DwMasterTr::terminationTest(int nnewcols, int itercnt, double relgap) {

	bool term = false;

	if (phase_ == 1) {
		return DwAlgo::terminationTest(nnewcols, itercnt, relgap);
	} else {
#if 0
		std::vector<int> delcols;
		for (int j = ncols_tr_, k = 0; j < si_->getNumCols(); ++k) {
			if (cols_generated_[k]->active_ == false)
				continue;
			if (fabs(si_->getColSolution()[j]) < 1.0e-8 && fabs(si_->getReducedCost()[j]) > 1.0e-8)
				cols_generated_[k]->age_++;
			else
				cols_generated_[k]->age_ = 0;
			if (cols_generated_[k]->age_ == 10) {
				cols_generated_[k]->active_ = false;
				delcols.push_back(j);
			}
			j++;
		}
		if (delcols.size() > 0) {
			CoinWarmStartBasis* ws = dynamic_cast<CoinWarmStartBasis*>(si_->getWarmStart());
			ws->deleteColumns(delcols.size(), &delcols[0]);
			si_->deleteCols(delcols.size(), &delcols[0]);
			si_->setWarmStart(ws);
			FREE_PTR(ws)
			message_->print(1, "  Deleted %u columns.\n", delcols.size());
		}
#endif
#if 0
		double tr_term = 0.0;
		double g = 0.0;
		for (int j = 0; j < nrows_orig_ + nrows_branch_; ++j) {
			tr_term += si_->getObjCoefficients()[2*j] * si_->getColSolution()[2*j];
			tr_term += si_->getObjCoefficients()[2*j+1] * si_->getColSolution()[2*j+1];
//		if (si_->getColSolution()[2*j] + si_->getColSolution()[2*j+1] > 1.0e-10)
//			DSPdebugMessage("row %d: obj [%e, %e], g [%e, %e]\n",
//					2*j, si_->getObjCoefficients()[2*j], si_->getObjCoefficients()[2*j+1],
//					si_->getColSolution()[2*j], si_->getColSolution()[2*j+1]);
			g += si_->getColSolution()[2*j] + si_->getColSolution()[2*j+1];
		}
		DSPdebugMessage("g %e, tr_term %e, master objective w/t tr_term %e\n", g, tr_term, primobj_ - tr_term);
#endif
		if (iterlim_ <= itercnt ||
			relgap < 0.01 - 1.0e-8 ||
			(phase_ == 1 && dualobj_ - bestprimobj_ > -1.0e-10)) {
			status_ = DSP_STAT_FEASIBLE;
			term = true;
#ifdef DSP_DEBUG
			for (unsigned k = 0; k < cols_generated_.size(); ++k)
				message_->print(1, "Column %d: age %d\n", k, cols_generated_[k]->age_);
#endif
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
			message_->print(2, "  [TR] SERIOUS STEP: dual objective %e", curLb);

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
			message_->print(2, "\n");
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

	double dist, mindist = 1.0;
	int find = -1; /**< fractional index */

	BGN_TRY_CATCH

	DSP_RTN_CHECK_RTN_CODE(
			preHeuristic(rlbd, rubd, primsol_org, primobj_org, dualobj_org, status_org));

#define DIVE_WITH_BACKTRACK
#ifdef DIVE_WITH_BACKTRACK
	std::vector<CoinTriple<int,int,double> > branchList;
	DSP_RTN_CHECK_RTN_CODE(gutsOfDive(branchList, 0));
#else
	status_ = DSP_STAT_FEASIBLE;
	while (status_ == DSP_STAT_FEASIBLE) {
		/** find a fractional */
		find = -1;
		mindist = 1.0;
		for (int i = 0; i < nrows_branch_; ++i) {
			int colind = branch_row_to_col_[nrows_orig_ + i];
			/** do not consider those with negative objective coefficients */
			if (org_obj_[colind] >= 0) continue;

			dist = fabs(primsol_[colind] - floor(primsol_[colind] + 0.5));
			if (dist < 1.0e-8) continue;
			if (dist < mindist) {
				mindist = dist;
				find = i;
			}
		}

		if (find > -1) {
			/** fix bound */
			double rounded = round(primsol_[branch_row_to_col_[nrows_orig_+find]]);
			si_->setRowBounds(nrows_orig_ + find, rounded, rounded);
			worker_->setColBounds(branch_row_to_col_[nrows_orig_ + find], rounded, rounded);
			DSPdebugMessage("Diving fixed variable %d to %e.\n", find, rounded);

			/** solve */
			DSP_RTN_CHECK_RTN_CODE(solvePhase1());

			if (primobj_ > 1.0e-8)
				status_ = DSP_STAT_PRIM_INFEASIBLE;
			else
				status_ = DSP_STAT_FEASIBLE;
		} else {

			/** find a fractional */
			for (int i = 0; i < nrows_branch_; ++i) {
				int colind = branch_row_to_col_[nrows_orig_ + i];
				/** do not consider those with negative objective coefficients */
				if (org_obj_[colind] < 0) continue;

				dist = fabs(primsol_[colind] - floor(primsol_[colind] + 0.5));
				if (dist < 1.0e-8) continue;
				if (dist < mindist) {
					mindist = dist;
					find = i;
				}
			}

			if (find > -1) {
				/** fix bound */
				double rounded = round(primsol_[branch_row_to_col_[nrows_orig_+find]]);
				si_->setRowBounds(nrows_orig_ + find, rounded, rounded);
				worker_->setColBounds(branch_row_to_col_[nrows_orig_ + find], rounded, rounded);
				DSPdebugMessage("Diving fixed variable %d to %e.\n", find, rounded);

				/** solve */
				DSP_RTN_CHECK_RTN_CODE(solvePhase1());

				if (primobj_ > 1.0e-8)
					status_ = DSP_STAT_PRIM_INFEASIBLE;
				else
					status_ = DSP_STAT_FEASIBLE;
			} else {
				DSPdebugMessage("Diving found no fractional variables.\n");
				DSP_RTN_CHECK_RTN_CODE(solvePhase2());
			}
		}

		/** collect solutions */
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
			if (find < 0) {
				/** check integer feasibility */
				bool fractional = false;
				for (int i = 0; i < nrows_branch_; ++i) {
					double x = primsol_[branch_row_to_col_[nrows_orig_ + i]];
					if (fabs(x - floor(x + 0.5)) > 1.0e-10) {
						fractional = true;
						DSPdebugMessage("Heuristic found a fractional solution (x %d [%e]).\n", i, x);
						break;
					}
				}
				if (!fractional) {
					DSPdebugMessage("Heuristic found an integer solution.\n");
					if (bestprimobj_ > primobj_) {
						bestprimobj_ = primobj_;
						CoinCopyN(primsol_, ncols_orig_, bestprimsol_);
						DSPdebugMessage("Heuristic updated the best upper bound %e.\n", bestprimobj_);
					}
					break;
				}
			}
		}
	}
#endif

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

	int status; /**< status at the current depth */

	BGN_TRY_CATCH

	if (depth > 0) {
		/** solve */
		DSP_RTN_CHECK_RTN_CODE(solvePhase1());

		/** determine solution status */
		if (primobj_ > 1.0e-8)
			status_ = DSP_STAT_PRIM_INFEASIBLE;
		else {
			status_ = DSP_STAT_FEASIBLE;

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
		DSPdebugMessage("Diving solution status %d\n", status_);
	}

	while (1) {
		if (status_ == DSP_STAT_PRIM_INFEASIBLE)
			break;

		/** find a fractional */
		findPhase = 0;
		branchIndex = -1;
		mindist = 1.0;
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

				double candValue = primsol_[colind];
				double candRounded = round(candValue);
				int candDirection = (candRounded < candValue) ? -1 : 1;

				/** check if the branch found is in branchList. */
				for (unsigned j = 0; j < branchList.size(); ++j) {
					if (branchList[j].first == i &&
						branchList[j].second == candDirection &&
						fabs(branchList[j].third - candRounded) < 1.0e-10) {
						/** yes, it exists; so do not branch on this. */
						dist = 1.0;
						break;
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
			DSPdebugMessage("Diving fixed variable %d [%e] to %e (depth %d).\n", branchIndex, branchValue, rounded, depth);

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
		} else {
			DSPdebugMessage("Diving found no fractional variables.\n");
			DSP_RTN_CHECK_RTN_CODE(solvePhase2());

			/** collect solutions */
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
				/** check integer feasibility */
				bool fractional = false;
				for (int i = 0; i < nrows_branch_; ++i) {
					double x = primsol_[branch_row_to_col_[nrows_orig_ + i]];
					if (fabs(x - floor(x + 0.5)) > 1.0e-10) {
						fractional = true;
						DSPdebugMessage("Heuristic found a fractional solution (x %d [%e]).\n", i, x);
						break;
					}
				}
				if (!fractional) {
					DSPdebugMessage("Heuristic found an integer solution.\n");
					primobj_ = 0.0;
					for (int j = ncols_tr_; j < si_->getNumCols(); ++j)
						primobj_ += si_->getObjCoefficients()[j] * si_->getColSolution()[j];
					if (bestprimobj_ > primobj_) {
						bestprimobj_ = primobj_;
						CoinCopyN(primsol_, ncols_orig_, bestprimsol_);
						DSPdebugMessage("Heuristic updated the best upper bound %e.\n", bestprimobj_);
					}
					break;
				}
			}
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

	message_->logLevel_ = 0;

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
