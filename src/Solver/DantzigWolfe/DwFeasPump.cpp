/*
 * DwFeasPump.cpp
 *
 *  Created on: Oct 27, 2016
 *      Author: kibaekkim
 */

#define DSP_DEBUG

#include <random>
/** Coin */
#include "OsiClpSolverInterface.hpp"
#include "OsiCpxSolverInterface.hpp"
/** Dsp */
#include "Model/TssModel.h"
#include "Solver/DantzigWolfe/DwFeasPump.h"
#include "Solver/DantzigWolfe/DwMaster.h"

DwFeasPump::DwFeasPump(DwMaster* master):
DwAlgo(master->getModelPtr(), master->getParPtr(), master->getMessagePtr()) {

	/** initialize random seed */
	srand(1);

	BGN_TRY_CATCH

	org_mat_ = new CoinPackedMatrix(*master->org_mat_);
	org_clbd_  = new double [org_mat_->getNumCols()];
	org_cubd_  = new double [org_mat_->getNumCols()];
	org_ctype_ = new char [org_mat_->getNumCols()];
	org_obj_   = new double [org_mat_->getNumCols()];
	org_rlbd_  = new double [org_mat_->getNumRows()];
	org_rubd_  = new double [org_mat_->getNumRows()];
	CoinCopyN(master->org_clbd_, org_mat_->getNumCols(), org_clbd_);
	CoinCopyN(master->org_cubd_, org_mat_->getNumCols(), org_cubd_);
	CoinCopyN(master->org_ctype_, org_mat_->getNumCols(), org_ctype_);
	CoinCopyN(master->org_obj_, org_mat_->getNumCols(), org_obj_);
	CoinCopyN(master->org_rlbd_, org_mat_->getNumRows(), org_rlbd_);
	CoinCopyN(master->org_rubd_, org_mat_->getNumRows(), org_rubd_);

	ncols_orig_ = org_mat_->getNumCols(); /**< number of columns in the original master */
	nrows_orig_ = org_mat_->getNumRows(); /**< number of rows in the original master */
	nrows_branch_ = 0; /**< number of branching rows in the restricted master */
	nrows_conv_ = model_->getNumSubproblems(); /**< number of convex combination rows in the restricted master */

	/** maps each branching row to original column index */
	for (int j = 0; j < ncols_orig_; ++j)
		if (org_ctype_[j] != 'C') {
			branch_row_to_col_[nrows_orig_ + nrows_branch_] = j;
			nrows_branch_++;
		}

	/** number of rows in the restricted master */
	nrows_ = nrows_orig_ + nrows_branch_ + nrows_conv_;
	DSPdebugMessage("nrwos_ %d, nrows_orig_ %d, nrows_branch_ %d, nrows_conv_ %d\n",
			nrows_, nrows_orig_, nrows_branch_, nrows_conv_);

	/** create subproblem solver */
    worker_ = new DwWorker(model_, par_, message_);

	/** create problem */
	DSP_RTN_CHECK(createProblem());

	END_TRY_CATCH(;)
}

DSP_RTN_CODE DwFeasPump::setBranchRowBounds(
		const double* rlbd,
		const double* rubd) {
	BGN_TRY_CATCH
	for (int i = 0; i < nrows_branch_; ++i) {
		si_->setRowBounds(nrows_orig_+i, rlbd[i], rubd[i]);
		rlbd_branch_[i] = rlbd[i];
		rubd_branch_[i] = rubd[i];
	}
	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

DSP_RTN_CODE DwFeasPump::createProblem() {
#define FREE_MEMORY      \
	FREE_PTR(mat)        \
	FREE_ARRAY_PTR(rlbd) \
	FREE_ARRAY_PTR(rubd)

	/** master problem */
	CoinPackedMatrix * mat = NULL;
	double * rlbd = NULL;
	double * rubd = NULL;

	BGN_TRY_CATCH

	/** allocate memory */
	rlbd = new double [nrows_];
	rubd = new double [nrows_];
	rlbd_branch_ = new double [nrows_branch_];
	rubd_branch_ = new double [nrows_branch_];

	/** create column-wise matrix and set number of rows */
	mat = new CoinPackedMatrix(true, 0, 0);
	mat->setDimensions(nrows_, 0);

	/** Set row bounds */
	CoinCopyN(org_rlbd_, nrows_orig_, rlbd);
	CoinCopyN(org_rubd_, nrows_orig_, rubd);
	for (int i = 0; i < nrows_branch_; ++i) {
		int j = branch_row_to_col_[nrows_orig_ + i];
		rlbd[nrows_orig_+i] = 0.0;
		rubd[nrows_orig_+i] = 0.0;
		rlbd_branch_[i] = org_clbd_[j];
		rubd_branch_[i] = org_cubd_[j];
	}
	CoinFillN(rlbd + nrows_orig_ + nrows_branch_, nrows_conv_, 1.0);
	CoinFillN(rubd + nrows_orig_ + nrows_branch_, nrows_conv_, 1.0);

	/** create solver */
	si_ = new OsiCpxSolverInterface();
	//dynamic_cast<OsiClpSolverInterface*>(si_)->getModelPtr()->setLogLevel(0);
	si_->messageHandler()->logLevel(0);
	//DSPdebug(si_->messageHandler()->logLevel(4));

	/** load problem data */
	si_->loadProblem(*mat, NULL, NULL, NULL, rlbd, rubd);
	DSPdebug(mat->verifyMtx(4));

	/** create auxiliary columns for Phase 1 problem */
	double auxcolval;
	for (int i = nrows_orig_; i < nrows_orig_ + nrows_branch_; ++i) {
		auxcolval = 1.0;
		auxcols_.push_back(new CoinPackedVector(1, &i, &auxcolval));
		auxcolval = -1.0;
		auxcols_.push_back(new CoinPackedVector(1, &i, &auxcolval));
	}
	phase_ = 2;

	/** allocate memory for solution */
	bestprimsol_ = new double [ncols_orig_];

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	/** release memory */
	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

DSP_RTN_CODE DwFeasPump::solve() {
#define FREE_MEMORY \
	FREE_ARRAY_PTR(activity) \
	FREE_PTR(rounded)

	/** original restricted mater problem information */
	double* activity = NULL;
	bool cycle = false;
	int cyclechk = 0;
	CoinPackedVector* rounded = NULL;

	int maxiter = 100;
	int maxcyclechk = 1000; /**< maximum number of checking cycles */
	int maxLookback = 1000000; /**< maximum number of iteration to look back for detecting a cycle */
	int maxFlips = 20;//floor(nrows_branch_ * 0.5); /**< maximum number of variables to flip when cycle is detected. */
	CoinPackedVector fracvals;

	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution(0.5,1.5);

	BGN_TRY_CATCH

	/** resolve  */
	si_->resolve();

	/** copy bounds */
	activity = new double [nrows_branch_];
	CoinCopyN(si_->getRowActivity() + nrows_orig_, nrows_branch_, activity);

	message_->print(1, "Running feasibility pump heuristic.\n");
	int iter = 1;
	while (iter < maxiter) {

		/** round and fix */
		rounded = new CoinPackedVector;
		for (int i = 0; i < nrows_branch_; ++i) {
			double rndval = round(activity[i]);
			si_->setRowBounds(nrows_orig_+i, rndval, rndval);
			if (fabs(rndval) > 1.0e-8)
				rounded->insert(i, rndval);
		}
		message_->print(2, "[FP] Iteration %3d, rounded and fixed fractional values.\n", iter);
		//DSPdebugMessage("Rounded integer value:\n");
		//DSPdebug(DspMessage::printArray(&rounded));

		/** detect cycle */
		int numVisited = visited_.size();
		if (numVisited > 0) {
			cyclechk = 0;
			cycle = true;
			while (cycle && cyclechk < maxcyclechk) {
				message_->print(3, "FP: Checking cycles.\n");
				cyclechk++;
				cycle = false;
				int v;
				for (v = numVisited - 1; v >= CoinMax(0,numVisited - maxLookback); --v) {
					if (rounded->getNumElements() == visited_[v]->getNumElements() &&
							std::equal(rounded->getIndices(), rounded->getIndices() + rounded->getNumElements(), visited_[v]->getIndices()) &&
							std::equal(rounded->getElements(), rounded->getElements() + rounded->getNumElements(), visited_[v]->getElements())) {
						message_->print(3, "FP: A cycle is detected.\n");
						cycle = true;
						break;
					}
					if (v == 0) break;
				}
				if (!cycle) break;

				/** if immediate cycle */
				if (v == numVisited - 1) {
					/** find and sort fractional variables */
					fracvals.clear();
					for (int i = 0; i < nrows_branch_; ++i) {
						double rndval = si_->getRowUpper()[nrows_orig_+i];
						double frac = fabs(rndval - activity[i]);
						fracvals.insert(i,frac);
					}
					fracvals.sortDecrElement();
					/** flip T most fractional variables */
					int maxFlipsRand = floor(maxFlips * distribution(generator));
					message_->print(3, "FP: Flip %d most fractional variables.\n", maxFlipsRand);
					for (int i = 0; i < fracvals.getNumElements(); ++i) {
						if (i >= maxFlipsRand) break;
						double rndval = si_->getRowUpper()[nrows_orig_+i];
						if (rndval == rubd_branch_[i] && rndval == rlbd_branch_[i])
							continue;
						else if (rndval == rubd_branch_[i])
							si_->setRowBounds(nrows_orig_+i, rndval-1, rndval-1);
						else if (rndval == rlbd_branch_[i])
							si_->setRowBounds(nrows_orig_+i, rndval+1, rndval+1);
						else if (distribution(generator) < 1.0)
							si_->setRowBounds(nrows_orig_+i, rndval-1, rndval-1);
						else
							si_->setRowBounds(nrows_orig_+i, rndval+1, rndval+1);
					}
				} else {
					/** randomize */
					message_->print(3, "FP: Randomly flip variables.\n");
					for (int i = 0; i < nrows_branch_; ++i) {
						double rndval = si_->getRowUpper()[nrows_orig_+i];
						if (rndval == rubd_branch_[i] && rndval == rlbd_branch_[i])
							continue;
						if (distribution(generator) < 1.0)
							continue;
						else if (rndval == rubd_branch_[i])
							si_->setRowBounds(nrows_orig_+i, rndval-1, rndval-1);
						else if (rndval == rlbd_branch_[i])
							si_->setRowBounds(nrows_orig_+i, rndval+1, rndval+1);
						else if (distribution(generator) < 1.0)
							si_->setRowBounds(nrows_orig_+i, rndval-1, rndval-1);
						else
							si_->setRowBounds(nrows_orig_+i, rndval+1, rndval+1);
					}
				}
				/** store solution */
				rounded->clear();
				for (int i = 0; i < nrows_branch_; ++i) {
					double rndval = si_->getRowLower()[nrows_orig_+i];
					if (fabs(rndval) > 1.0e-8)
						rounded->insert(i, rndval);
				}
			}
		}

		/** reached maximum cycle checks */
		if (cycle) {
			message_->print(2, "FP: Maximum cycle checks %d reached.\n", cyclechk);
			break;
		}

		visited_.push_back(rounded);
		rounded = NULL;
		message_->print(3, "FP: %u solutions visited.\n", visited_.size());

		/** resolve FP master */
		si_->resolve();

		/** TODO:
		 * When to generate columns??? We do not know if the problem will end up being infeasible a priori,
		 * whereas the FP problem is always feasible.
		 */

		/** solve projection LP */
		DSP_RTN_CHECK_RTN_CODE(solvePhase1());

		if (si_->getObjValue() < 1.0e-8) {

			/** Improve upper bound */
			DSP_RTN_CHECK_RTN_CODE(solvePhase2());

			/** get upper bound */
			double ub = si_->getObjValue();
			message_->print(2, "FP: Found a primal bound %e.\n", ub);

			if (ub < bestprimobj_) {
				message_->print(2, "FP: Update best primal bound %e.\n", ub);
				bestprimobj_ = ub;

				/** recover original solution */
				CoinZeroN(bestprimsol_, ncols_orig_);
				for (unsigned k = 0, j = 0; k < cols_generated_.size(); ++k)
					/** consider active columns only */
					if (cols_generated_[k]->active_) {
						CoinPackedVector xlam = cols_generated_[k]->x_ * si_->getColSolution()[j];
						for (int i = 0; i < xlam.getNumElements(); ++i)
							bestprimsol_[xlam.getIndices()[i]] += xlam.getElements()[i];
						j++;
					}
			}
			break;
		}
		message_->print(3, "FP: Projected to the LP relaxation (infeasibility %e).\n", si_->getObjValue());

		/** copy integer variable value */
		CoinCopyN(si_->getRowActivity() + nrows_orig_, nrows_branch_, activity);
		const CoinPackedMatrix* matByCol = si_->getMatrixByCol();
		for (unsigned j = 0; j < auxcolindices_.size(); ++j) {
			for (int i = 0; i < nrows_branch_; ++i)
				activity[i] -= matByCol->getCoefficient(nrows_orig_+i, auxcolindices_[j]) * si_->getColSolution()[auxcolindices_[j]];
		}

		/** increment iteration count */
		iter++;
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

DSP_RTN_CODE DwFeasPump::copyColumns(std::vector<DwCol*>& cols) {
	BGN_TRY_CATCH

	/** remove all columns */
	std::vector<int> delcols;
	for (int j = 0; j < si_->getNumCols(); ++j)
		delcols.push_back(j);
	si_->deleteCols(delcols.size(), &delcols[0]);
	auxcolindices_.clear();
	phase_ = 2;

	/** add active columns */
	for (unsigned k = 0; k < cols.size(); ++k) {
		if (cols[k]->active_)
			si_->addCol(cols[k]->col_, cols[k]->lb_, cols[k]->ub_, 0.0);
		cols_generated_.push_back(cols[k]);
	}
	DSPdebugMessage("Master has %u columns generated. FP has %u columns copied.\n",
			cols.size(), cols_generated_.size());

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DwFeasPump::getNewCols(std::vector<DwCol*>& cols) {
	BGN_TRY_CATCH

	/** add active columns */
	DSPdebugMessage("Master has %u columns generated. FP has %u columns generated.\n",
			cols.size(), cols_generated_.size());
	for (unsigned k = 0; k < cols_generated_.size(); ++k) {
		if (k < cols.size()) {
			cols_generated_[k] = NULL;
			continue;
		}
		cols.push_back(cols_generated_[k]);
		cols_generated_[k] = NULL;
	}
	cols_generated_.clear();

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}
