/*
 * DwSolver.cpp
 *
 *  Created on: Oct 27, 2016
 *      Author: kibaekkim
 */

//#define DSP_DEBUG

/** Coin */
#include "OsiCpxSolverInterface.hpp"
/** Dsp */
#include "Utility/DspUtility.h"
#include "Solver/DantzigWolfe/DwAlgo.h"
#include "Model/TssModel.h"

DwAlgo::DwAlgo(
		DecModel *   model,  /**< model pointer */
        DspParams *  par,    /**< parameters */
        DspMessage * message /**< message pointer */):
DecSolver(model, par, message),
useCpxBarrier_(false),
phase_(1),
ncols_orig_(0),
nrows_(0),
nrows_orig_(0),
nrows_branch_(0),
nrows_conv_(0),
rlbd_branch_(NULL),
rubd_branch_(NULL),
worker_(NULL),
org_mat_(NULL),
org_clbd_(NULL),
org_cubd_(NULL),
org_obj_(NULL),
org_ctype_(NULL),
org_rlbd_(NULL),
org_rubd_(NULL) {
}

DwAlgo::~DwAlgo() {
	for (unsigned i = 0; i < auxcols_.size(); ++i)
		FREE_PTR(auxcols_[i]);
	FREE_ARRAY_PTR(rlbd_branch_);
	FREE_ARRAY_PTR(rubd_branch_);
	FREE_PTR(worker_);
	FREE_PTR(org_mat_);
	FREE_ARRAY_PTR(org_clbd_);
	FREE_ARRAY_PTR(org_cubd_);
	FREE_ARRAY_PTR(org_obj_);
	FREE_ARRAY_PTR(org_ctype_);
	FREE_ARRAY_PTR(org_rlbd_);
	FREE_ARRAY_PTR(org_rubd_);
	for (unsigned i = 0; i < cols_generated_.size(); ++i)
		FREE_PTR(cols_generated_[i]);
}

DSP_RTN_CODE DwAlgo::solve() {
	BGN_TRY_CATCH

	if (phase_ == 2) {
		DSP_RTN_CHECK_RTN_CODE(solvePhase2());
		if (status_ == DSP_STAT_PRIM_INFEASIBLE) {
			DSPdebugMessage("Converting to Phase 1.\n");
			DSP_RTN_CHECK_RTN_CODE(solvePhase1());
		}
	} else
		DSP_RTN_CHECK_RTN_CODE(solvePhase1());

	if (phase_ == 1) {
		if (status_ == DSP_STAT_FEASIBLE || status_ == DSP_STAT_OPTIMAL) {
			if (primobj_ > 1.0e-8)
				status_ = DSP_STAT_PRIM_INFEASIBLE;
			else {
				DSPdebugMessage("Converting to Phase 2.\n");
				DSP_RTN_CHECK_RTN_CODE(solvePhase2());
			}
		}
	}

	/** collect solutions */
	if (status_ == DSP_STAT_OPTIMAL || status_ == DSP_STAT_FEASIBLE) {
#if 1
		/** recover original solution */
		CoinZeroN(primsol_, ncols_orig_);
		for (unsigned k = 0, j = 0; k < cols_generated_.size(); ++k) {
			/** do not consider inactive columns */
			if (cols_generated_[k]->active_ == false)
				continue;
			CoinPackedVector xlam = cols_generated_[k]->x_ * si_->getColSolution()[j];
			for (int i = 0; i < xlam.getNumElements(); ++i)
				primsol_[xlam.getIndices()[i]] += xlam.getElements()[i];
			j++;
		}
#endif

		/** heuristics */
		DSP_RTN_CHECK_RTN_CODE(heuristics());
	} else if (status_ == DSP_STAT_DUAL_INFEASIBLE) {
		//DSPdebug(si_->writeMps("master"));
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DwAlgo::initialColumns() {
#define FREE_MEMORY     \
	FREE_ARRAY_PTR(piA) \
	FREE_ARRAY_PTR(Ax)

	/** column generation info */
	double* piA = NULL;
	double* Ax = NULL;
	std::vector<int> subinds;
	std::vector<int> substatus;
	std::vector<double> subcxs;
	std::vector<double> subobjs;
	std::vector<CoinPackedVector*> subsols;

	BGN_TRY_CATCH

	/** allocate memory */
	piA = new double [ncols_orig_];
	Ax  = new double [nrows_orig_];

	/** generate columns */
	CoinZeroN(piA, ncols_orig_);
	worker_->generateCols(2, piA, subinds, substatus, subcxs, subobjs, subsols);
	dualobj_ = 0.0;

	for (unsigned int s = 0; s < subinds.size(); ++s) {

		int sind = subinds[s];                  /**< actual subproblem index */
		const CoinPackedVector* x = subsols[s]; /**< retrieve subproblem solution */
		double newcoef = subobjs[s];            /**< column objective coefficient */
		//DSPdebugMessage("subproblem %d objective value %e\n", sind, newcoef);

		/** calculate initial dual bound */
		dualobj_ += subobjs[s];

		/** take A x^k */
		org_mat_->times(*x, Ax);

		/** create a column vector */
		CoinPackedVector colvec;
		colvec.reserve(nrows_);

		/** original constraints */
		for (int i = 0; i < nrows_orig_; ++i) {
			if (fabs(Ax[i]) > 1.0e-10)
				colvec.insert(i, Ax[i]);
		}

		/** branching constraints */
		for (int i = 0; i < nrows_branch_; ++i) {
			int j = branch_row_to_col_[nrows_orig_+i];
			int sparse_index = x->findIndex(j);
			if (sparse_index == -1) continue;
			double val = x->getElements()[sparse_index];
			if (fabs(val) > 1.0e-10)
				colvec.insert(nrows_orig_+i, val);
		}

		/** convex combination constraints, unless unbounded */
		if (substatus[s] != DSP_STAT_DUAL_INFEASIBLE)
			colvec.insert(nrows_orig_ + nrows_branch_ + sind, 1.0);
		//DSPdebugMessage("new column from subproblem %d:\n", sind);
		//DSPdebug(DspMessage::printArray(&colvec));

		/** store columns */
		CoinPackedVector copy_of_x(*x);
		cols_generated_.push_back(new DwCol(sind, copy_of_x, colvec, newcoef, 0.0, COIN_DBL_MAX));
	}
	message_->print(0, "Generated %u initial columns. Initial dual bound %e\n", cols_generated_.size(), dualobj_);

	/** clear column generation information */
	for (unsigned i = 0; i < subsols.size(); ++i)
		FREE_PTR(subsols[i]);

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

DSP_RTN_CODE DwAlgo::solvePhase1() {
	BGN_TRY_CATCH

	if (phase_ == 2) {
		/** set phase 1 problem */
		int ncols = si_->getNumCols();
		for (int j = 0; j < ncols; ++j)
			si_->setObjCoeff(j, 0.0);
		for (unsigned j = 0; j < auxcols_.size(); ++j) {
			si_->addCol(*auxcols_[j], 0.0, COIN_DBL_MAX, 1.0);
			auxcolindices_.push_back(ncols + j);
		}
		phase_ = 1;
		DSPdebugMessage("Phase 1 has %d rows and %d columns.\n", si_->getNumRows(), si_->getNumCols());
	}

	DSP_RTN_CHECK_RTN_CODE(gutsOfSolve());

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

DSP_RTN_CODE DwAlgo::solvePhase2() {
	BGN_TRY_CATCH

	if (phase_ == 1) {
		/** set phase 2 problem */
		//DSPdebug(si_->writeMps("beforePhase2"));
		si_->deleteCols(auxcolindices_.size(), &auxcolindices_[0]);
		auxcolindices_.clear();
		DSPdebugMessage("Phase 2 has %d rows and %d columns.\n", si_->getNumRows(), si_->getNumCols());
		for (unsigned k = 0, j = 0; k < cols_generated_.size(); ++k)
			if (cols_generated_[k]->active_) {
				if (j >= si_->getNumCols()) {
					message_->print(0, "Trying to access invalid column index %d (ncols %d)\n", j, si_->getNumCols());
					return DSP_RTN_ERR;
				}
				si_->setObjCoeff(j, cols_generated_[k]->obj_);
				j++;
			}
		phase_ = 2;
		//DSPdebug(si_->writeMps("initialPhase2"));
	}

	DSP_RTN_CHECK_RTN_CODE(gutsOfSolve());
	//DSPdebug(si_->writeMps("afterPhase2"));

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

DSP_RTN_CODE DwAlgo::gutsOfSolve() {
#define FREE_MEMORY       \
	FREE_PTR(ws)          \
	FREE_ARRAY_PTR(price) \
	FREE_ARRAY_PTR(piA)

	OsiCpxSolverInterface* cpx = NULL;
	CoinWarmStartBasis* ws = NULL;
	double* price = NULL;
	double* piA = NULL;
	bool hasBasis = false;

	double feastol = 1.0e-8; /**< feasibility tolerance */
	double curLb = -COIN_DBL_MAX;

	/** column generation info */
	std::vector<int> subinds;
	std::vector<int> substatus;
	std::vector<double> subcxs;
	std::vector<double> subobjs;
	std::vector<CoinPackedVector*> subsols;

	BGN_TRY_CATCH

	cpx = dynamic_cast<OsiCpxSolverInterface*>(si_);

	/** allocate memory */
	price = new double [si_->getNumRows()];
	piA = new double [ncols_orig_];

	std::vector<double> prevsol(si_->getNumCols(), 0.0);
	//dualobj_ = -COIN_DBL_MAX;

	/** timing results */
	double stime;
	double t_total = CoinGetTimeOfDay();
	double t_master = 0.0;
	double t_gencols = 0.0;

	/** use dual simplex after branching */
	if (!useCpxBarrier_ || cpx == NULL) {
		si_->setHintParam(OsiDoDualInResolve, true);
		si_->messageHandler()->setLogLevel(0);
		hasBasis = true;
	}

	/** resolve */
	stime = CoinGetTimeOfDay();
	si_->resolve();
	primobj_ = si_->getObjValue();
	t_master += CoinGetTimeOfDay() - stime;
	convertCoinToDspStatus(si_, status_);

	/** get price */
	CoinCopyN(si_->getRowPrice(), si_->getNumRows(), price);

	/** get warm-start information */
	if (hasBasis)
		DSP_RTN_CHECK_RTN_CODE(getWarmStartInfo(prevsol, ws));

	/** use primal simplex after column generation */
	if (!useCpxBarrier_)
		si_->setHintParam(OsiDoDualInResolve, false);

	int itercnt = 0;

	while (status_ == DSP_STAT_OPTIMAL) {

		/** calculate pi^T A */
		DSP_RTN_CHECK_RTN_CODE(calculatePiA(price, piA));

		/** generate columns */
		stime = CoinGetTimeOfDay();
		DSP_RTN_CHECK_RTN_CODE(
				worker_->generateCols(phase_, piA, subinds, substatus, subcxs, subobjs, subsols));
		t_gencols += CoinGetTimeOfDay() - stime;

		/** termination test */
		if (terminationTestColgen(substatus))
			break;

		/** calculate lower bound */
		if (phase_ == 2) {
			DSP_RTN_CHECK_RTN_CODE(getLagrangianBound(price, subobjs, curLb));
			DSPdebugMessage("Current lower bound %e, best lower bound %e\n", curLb, dualobj_);
		}

		/** create and add columns */
		int nColsAdded = 0;
		DSP_RTN_CHECK_RTN_CODE(
				addCols(price, piA, subinds, substatus, subcxs, subobjs, subsols, nColsAdded));

		/** update master */
		DSP_RTN_CHECK_RTN_CODE(updateModel(price, curLb));

		/** set warm-start information */
		if (hasBasis)
			DSP_RTN_CHECK_RTN_CODE(setWarmStartInfo(prevsol, ws));

		/** re-optimize the master */
		stime = CoinGetTimeOfDay();
		si_->resolve();
		primobj_ = si_->getObjValue();
		t_master += CoinGetTimeOfDay() - stime;
		convertCoinToDspStatus(si_, status_);

		/** get price */
		CoinCopyN(si_->getRowPrice(), si_->getNumRows(), price);

		double relgap = (primobj_-dualobj_)/(1.0e-10+fabs(primobj_))*100;
		message_->print(0, "[Phase %d] Iteration %3d: Master objective %e, ", phase_, itercnt, primobj_);
		if (phase_ == 2)
			message_->print(0, "Lb %e (gap %.2f %%), ", dualobj_, relgap);
		message_->print(0, "nrows %d, ncols %d (new cols %d), ", si_->getNumRows(), si_->getNumCols(), nColsAdded);
		message_->print(0, "itercnt %d, timing (total %.2f, master %.2f, gencols %.2f), statue %d\n",
				si_->getIterationCount(), CoinGetTimeOfDay() - t_total, t_master, t_gencols, status_);

		/** termination test */
		if (terminationTest(nColsAdded, itercnt, relgap))
			break;

		/** get warm-start information */
		if (hasBasis)
			DSP_RTN_CHECK_RTN_CODE(getWarmStartInfo(prevsol, ws));

		itercnt++;

#ifdef DSP_DEBUG1
		char fname[128];
		sprintf(fname, "master%d", itercnt);
		si_->writeMps(fname);
#endif
	}

	/** free memory for subproblem solutions */
	for (unsigned i = 0; i < subsols.size(); ++i)
		FREE_PTR(subsols[i]);
	subsols.clear();

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

DSP_RTN_CODE DwAlgo::addCols(
		const double* price,                  /**< [in] price */
		const double* piA,                    /**< [in] pi^T A */
		std::vector<int>& indices,            /**< [in] subproblem indices corresponding to cols*/
		std::vector<int>& statuses,           /**< [in] subproblem solution status */
		std::vector<double>& cxs,             /**< [in] solution times original objective coefficients */
		std::vector<double>& objs,            /**< [in] subproblem objective values */
		std::vector<CoinPackedVector*>& sols, /**< [in] subproblem solutions */
		int& nadded                           /**< [out] number of columns added */) {
#define FREE_MEMORY \
	FREE_ARRAY_PTR(Ax)

	double* Ax = NULL;
	CoinPackedVector colvec;
	TssModel* tss = NULL;

	BGN_TRY_CATCH

	/** allocate memory */
	Ax = new double [nrows_orig_];
	colvec.reserve(nrows_);

	if (model_->isStochastic())
		tss = dynamic_cast<TssModel*>(model_);

	/** reset counter */
	nadded = 0;

	for (unsigned int s = 0; s < indices.size(); ++s) {
		int sind = indices[s]; /**< actual subproblem index */

		/** cutoff = dual variable corresponding to the convex-combination constraint */
		double cutoff = price[nrows_orig_ + nrows_branch_ + sind];
		DSPdebugMessage("pricing out: %e < %e ? (colobj %e, status %d)\n", objs[s], cutoff, cxs[s], statuses[s]);

		if (statuses[s] == DSP_STAT_DUAL_INFEASIBLE || objs[s] < cutoff - 1.0e-4) {
			/** retrieve subproblem solution */
			const CoinPackedVector* x = sols[s];

			/** create a column objective */
			double newcoef = cxs[s];

			/** take A x^k */
			org_mat_->times(*x, Ax);

			/** clear a column vector */
			colvec.clear();

			/** original constraints */
			for (int i = 0; i < nrows_orig_; ++i)
				if (fabs(Ax[i]) > 1.0e-10)
					colvec.insert(i, Ax[i]);

			/** skip zero column */
//			if (colvec.getNumElements() == 0)
//				continue;

			/** branching constraints */
			for (int i = 0; i < nrows_branch_; ++i) {
				int j = branch_row_to_col_[nrows_orig_+i];
				int sparse_index = x->findIndex(j);
				if (sparse_index == -1) continue;
				double val = x->getElements()[sparse_index];
				if (fabs(val) > 1.0e-10)
					colvec.insert(nrows_orig_ + i, val);
			}

			/** convex combination constraints */
			if (statuses[s] != DSP_STAT_DUAL_INFEASIBLE)
				colvec.insert(nrows_orig_ + nrows_branch_ + sind, 1.0);
//			DSPdebugMessage("new column from subproblem %d (objcoef %e):\n", sind, newcoef);
//			DSPdebug(DspMessage::printArray(&colvec));

			/** add the column vector */
			if (phase_ == 1)
				si_->addCol(colvec, 0.0, COIN_DBL_MAX, 0.0);
			else if (phase_ == 2)
				si_->addCol(colvec, 0.0, COIN_DBL_MAX, newcoef);

			/** store columns */
			cols_generated_.push_back(new DwCol(sind, *x, colvec, newcoef, 0.0, COIN_DBL_MAX));
			nadded++;
		}
	}
	//message_->print(3, "Number of columns in the pool: %u\n", cols_generated_.size());

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

DSP_RTN_CODE DwAlgo::getLagrangianBound(
		const double* price,       /**< [in] price */
		std::vector<double>& objs, /**< [in] subproblem objective values */
		double& lb                 /**< [out] lower bound */) {
	BGN_TRY_CATCH

	/** calculate lower bound */
	lb = 0;
	for (int j = 0; j < nrows_orig_ + nrows_branch_; ++j) {
		if (fabs(price[j]) < 1.0e-8) continue;
		switch (si_->getRowSense()[j]) {
		case 'E':
			lb += price[j] * si_->getRowUpper()[j];
			break;
		case 'L':
			lb += price[j] * si_->getRowUpper()[j];
			break;
		case 'G':
			lb += price[j] * si_->getRowLower()[j];
			break;
		case 'R':
			if (price[j] > 0)
				lb += price[j] * si_->getRowLower()[j];
			else
				lb += price[j] * si_->getRowUpper()[j];
			break;
		}
	}

	for (unsigned int s = 0; s < objs.size(); ++s)
		lb += objs[s];

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

bool DwAlgo::terminationTest(int nnewcols, int itercnt, double relgap) {

	if (nnewcols == 0)
		true;

	bool term = false;
	double feastol = 1.0e-8;

	if (phase_ == 1) {
		if (si_->isProvenOptimal() && primobj_ < feastol) {
			status_ = DSP_STAT_FEASIBLE;
			DSPdebugMessage("Phase 1 found a feasible solution! %e < %e\n", si_->getObjValue(), feastol);
			term = true;
		}
	}
	if (phase_ == 2) {
		if (iterlim_ <= itercnt ||
			relgap < 0.01 ||
			dualobj_ >= bestprimobj_) {
			status_ = DSP_STAT_FEASIBLE;
			term = true;
		}
	}

	return term;
}

DSP_RTN_CODE DwAlgo::updateModel(
		const double* price, /**< [in] price */
		double curLb         /**< [in] current lower bound */) {
	BGN_TRY_CATCH

	/** update the best dual objective */
	if (dualobj_ < curLb) dualobj_ = curLb;

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

DSP_RTN_CODE DwAlgo::calculatePiA(
		const double* price, /**< [in] price */
		double*& piA         /**< [out] pi^T A */) {
	BGN_TRY_CATCH

	org_mat_->transposeTimes(price, piA);
	for (int i = nrows_orig_; i < nrows_orig_ + nrows_branch_; ++i)
		piA[branch_row_to_col_[i]] += price[i];
//	DSPdebugMessage("piA (%d):\n", ncols_orig_);
//	DSPdebug(DspMessage::printArray(ncols_orig_, piA));

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

bool DwAlgo::terminationTestColgen(std::vector<int>& statuses) {
	bool term = false;
	//if (phase_ == 1) {
		for (unsigned s = 0; s < statuses.size(); ++s)
			if (statuses[s] != DSP_STAT_OPTIMAL && statuses[s] != DSP_STAT_DUAL_INFEASIBLE) {
				status_ = DSP_STAT_PRIM_INFEASIBLE;
				break;
			}
		if (status_ == DSP_STAT_PRIM_INFEASIBLE) {
			DSPdebugMessage("Subproblems are infeasible.\n");
			term = true;
		}
	//}
	return term;
}

DSP_RTN_CODE DwAlgo::getWarmStartInfo(
		std::vector<double>& sol, /**< [out] current solution */
		CoinWarmStartBasis*& ws   /**< [out] warmstart basis */) {
	BGN_TRY_CATCH

	/** store previous solution and basis */
	CoinCopyN(si_->getColSolution(), si_->getNumCols(), &sol[0]);
	ws = dynamic_cast<CoinWarmStartBasis*>(si_->getWarmStart());

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DwAlgo::setWarmStartInfo(
		std::vector<double>& sol, /**< [out] current solution */
		CoinWarmStartBasis*& ws   /**< [out] warmstart basis */) {
	BGN_TRY_CATCH

	/** set previous solution and basis */
	sol.resize(si_->getNumCols(), 0.0);
	si_->setColSolution(&sol[0]);

	ws->resize(si_->getNumRows(), si_->getNumCols());
	si_->setWarmStart(ws);
	FREE_PTR(ws);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}
