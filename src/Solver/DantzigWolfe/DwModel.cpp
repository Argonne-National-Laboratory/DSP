/*
 * DwModel.cpp
 *
 *  Created on: Feb 13, 2017
 *      Author: kibaekkim
 */

//#define DSP_DEBUG

#include <DantzigWolfe/DwModel.h>
#include <DantzigWolfe/DwHeuristic.h>
#include "Model/TssModel.h"

DwModel::DwModel(): DspModel(), master_(NULL), infeasibility_(0.0) {}

DwModel::DwModel(DecSolver* solver): DspModel(solver), infeasibility_(0.0) {
	master_ = dynamic_cast<DwMaster*>(solver_);
	primsol_.resize(master_->ncols_orig_);
	//heuristics_.push_back(new DwRounding("Rounding", *this));
	heuristics_.push_back(new DwSmip("Smip", *this));
}

DwModel::~DwModel() {
	for (unsigned i = 0; i < heuristics_.size(); ++i)
		delete heuristics_[i];
}

DSP_RTN_CODE DwModel::solve() {
	BGN_TRY_CATCH

	/** set best primal objective value */
	solver_->setBestPrimalObjective(bestprimobj_);

	/** solve master */
	solver_->solve();

	status_ = solver_->getStatus();

	switch (status_) {
	case DSP_STAT_OPTIMAL:
	case DSP_STAT_FEASIBLE:
	case DSP_STAT_LIM_ITERorTIME: {

		primobj_ = master_->getPrimalObjective();
		dualobj_ = -master_->getBestDualObjective();

		if (primobj_ < 1.0e+20) {
			/** parse solution */
			int cpos = 0;
			std::fill(primsol_.begin(), primsol_.begin() + master_->ncols_orig_, 0.0);
			for (auto it = master_->cols_generated_.begin(); it != master_->cols_generated_.end(); it++) {
				if ((*it)->active_) {
					if (fabs(master_->getPrimalSolution()[cpos]) > 1.0e-10) {
						for (int i = 0; i < (*it)->x_.getNumElements(); ++i) {
							if ((*it)->x_.getIndices()[i] < master_->ncols_orig_)
								primsol_[(*it)->x_.getIndices()[i]] += (*it)->x_.getElements()[i] * master_->getPrimalSolution()[cpos];
						}
					}
					cpos++;
				}
			}
			//DspMessage::printArray(cpos, master_->getPrimalSolution());

			/** calculate infeasibility */
			infeasibility_ = 0.0;
			for (int j = 0; j < master_->ncols_orig_; ++j)
				if (master_->ctype_orig_[j] != 'C') {
					infeasibility_ += fabs(primsol_[j] - floor(primsol_[j] + 0.5));
				}
			solver_->getMessagePtr()->print(3, "Infeasibility: %+e\n", infeasibility_);

			bool isViolated = false;
			for (int j = 0; j < master_->ncols_orig_; ++j) {
				double viol = std::max(master_->clbd_node_[j] - primsol_[j], primsol_[j] - master_->cubd_node_[j]);
				if (viol > 1.0e-6) {
					printf("Violated variable at %d by %e (%+e <= %+e <= %+e)\n", j, viol,
							master_->clbd_node_[j], primsol_[j], master_->cubd_node_[j]);
					isViolated = true;
				}
			}
			if (isViolated) throw "Invalid branching was performed.";

			/** run heuristics */
			if (par_->getBoolParam("DW/HEURISTICS") && infeasibility_ > 1.0e-6) {
				/** FIXME */
				bestprimobj_ = COIN_DBL_MAX;
				for (auto it = heuristics_.begin(); it != heuristics_.end(); it++) {
					solver_->getMessagePtr()->print(1, "Running [%s] heuristic:\n", (*it)->name());
					int found = (*it)->solution(bestprimobj_, bestprimsol_);
					//printf("found %d bestprimobj %+e\n", found, bestprimobj_);
				}
			}
		}

		break;
	}
	default:
		break;
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

bool DwModel::chooseBranchingObjects(
		DspBranch*& branchingUp, /**< [out] branching-up object */
		DspBranch*& branchingDn  /**< [out] branching-down object */) {
	int findPhase = 0;
	bool branched = false;
	double dist, maxdist = 1.0e-6;
	int branchingIndex = -1;
	double branchingValue;

	/** smip branching */
	int ncols_first_stage = -1;   /**< number of first-stage columns in dd form */
	int branchingFirstStage = -1; /**< branching index in first stage */
	TssModel* tss = NULL;

	BGN_TRY_CATCH

	/** cleanup */
	FREE_PTR(branchingUp)
	FREE_PTR(branchingDn)

	if (solver_->getModelPtr()->isStochastic()) {
		/** two-stage stochastic model */
		tss = dynamic_cast<TssModel*>(solver_->getModelPtr());
		ncols_first_stage = tss->getNumScenarios() * tss->getNumCols(0);
	}
	DSPdebugMessage("ncols_first_stage %d\n", ncols_first_stage);

	findPhase = 0;
	while (findPhase < 2 && branchingIndex < 0) {
		/** most fractional value */
		for (int j = 0; j < master_->ncols_orig_; ++j) {
			if (findPhase == 0 && j > ncols_first_stage)
				break;
			if (master_->ctype_orig_[j] == 'C') continue;
			dist = fabs(primsol_[j] - floor(primsol_[j] + 0.5));
			if (dist > maxdist) {
				maxdist = dist;
				branchingIndex = j;
				branchingValue = primsol_[j];
			}
		}

		/** for the first pass of smip, look through expected first-stage integer variable values */
		if (ncols_first_stage > 0 && findPhase == 0 && branchingIndex < 0) {
			maxdist = 1.0e-8;
			for (int j = 0; j < tss->getNumCols(0); ++j) {
				if (master_->ctype_orig_[j] == 'C') continue;
				double expval = 0.0;
				for (int s = 0; s < tss->getNumScenarios(); ++s)
					expval += primsol_[tss->getNumCols(0) * s + j] / tss->getNumScenarios();
				dist = fabs(expval - floor(expval + 0.5));
				if (dist > maxdist) {
					maxdist = dist;
					branchingIndex = j;
					branchingValue = primsol_[j];
				}
			}
		}

		findPhase++;
	}

	if (branchingIndex > -1) {

		branched = true;

		/** get branching index in first stage */
		if (branchingIndex < ncols_first_stage)
			branchingFirstStage = branchingIndex % tss->getNumCols(0);

		/** creating branching objects */
		branchingUp = new DspBranch();
		branchingDn = new DspBranch();
		for (int j = 0; j < master_->ncols_orig_; ++j) {
			if (master_->ctype_orig_[j] == 'C') continue;
			/** NOTE: branching on all the first-stage variables if SMIP */
			if (branchingIndex == j || (tss != NULL && branchingFirstStage == j % tss->getNumCols(0))) {
				DSPdebugMessage("Creating branch objects on column %d (value %e): [%e,%e] and [%e,%e]\n", 
					j, branchingValue, ceil(branchingValue), master_->cubd_node_[j], master_->clbd_node_[j], floor(branchingValue));
				branchingUp->push_back(j, ceil(branchingValue), master_->cubd_node_[j]);
				branchingDn->push_back(j, master_->clbd_node_[j], floor(branchingValue));
			} else if (master_->clbd_node_[j] > master_->clbd_orig_[j] || master_->cubd_node_[j] < master_->cubd_orig_[j]) {
				/** store any bound changes made in parent nodes */
				branchingUp->push_back(j, master_->clbd_node_[j], master_->cubd_node_[j]);
				branchingDn->push_back(j, master_->clbd_node_[j], master_->cubd_node_[j]);
			}
		}
		branchingUp->bestBound_ = master_->getBestDualObjective();
		branchingDn->bestBound_ = master_->getBestDualObjective();
		branchingUp->dualsol_.assign(master_->getBestDualSolution(), master_->getBestDualSolution() + master_->nrows_);
		branchingDn->dualsol_.assign(master_->getBestDualSolution(), master_->getBestDualSolution() + master_->nrows_);
	} else {
		DSPdebugMessage("No branch object is found.\n");
	}

	END_TRY_CATCH_RTN(;,false)

	return branched;
}

