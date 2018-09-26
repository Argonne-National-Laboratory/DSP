/*
 * DwBranchInt.cpp
 *
 *  Created on: May 1, 2018
 *      Author: Kibaek Kim
 */

//#define DSP_DEBUG

#include <DantzigWolfe/DwBranchInt.h>
#include <Model/TssModel.h>

bool DwBranchInt::chooseBranchingObjects(
		std::vector<DspBranchObj*>& branchingObjs /**< [out] branching objects */) {
	int findPhase = 0;
	bool branched = false;
	double dist, maxdist = 1.0e-6;
	int branchingIndex = -1;
	double branchingValue;

	DecSolver* solver = model_->getSolver();
	DwMaster* master = dynamic_cast<DwMaster*>(solver);
	std::vector<double> primsol = model_->getPrimalSolution();

	/** smip branching */
	int ncols_first_stage = -1;   /**< number of first-stage columns in dd form */
	int branchingFirstStage = -1; /**< branching index in first stage */
	TssModel* tss = NULL;

	DspBranchObj* branchingUp = NULL;
	DspBranchObj* branchingDn = NULL;

	BGN_TRY_CATCH

	if (solver->getModelPtr()->isStochastic()) {
		/** two-stage stochastic model */
		tss = dynamic_cast<TssModel*>(solver->getModelPtr());
		ncols_first_stage = tss->getNumScenarios() * tss->getNumCols(0);
	}
	DSPdebugMessage("ncols_first_stage %d\n", ncols_first_stage);

	findPhase = 0;
	while (findPhase < 2 && branchingIndex < 0) {
		/** most fractional value */
		for (int j = 0; j < master->ncols_orig_; ++j) {
			if (findPhase == 0 && j > ncols_first_stage)
				break;
			if (master->ctype_orig_[j] == 'C') continue;
			dist = fabs(primsol[j] - floor(primsol[j] + 0.5));
			if (dist > maxdist) {
				maxdist = dist;
				branchingIndex = j;
				branchingValue = primsol[j];
			}
		}

#if 0
		/** for the first pass of smip, look through expected first-stage integer variable values */
		if (ncols_first_stage > 0 && findPhase == 0 && branchingIndex < 0) {
			maxdist = 1.0e-6;
			for (int j = 0; j < tss->getNumCols(0); ++j) {
				if (master->ctype_orig_[j] == 'C') continue;
				double expval = 0.0;
				for (int s = 0; s < tss->getNumScenarios(); ++s)
					expval += primsol[tss->getNumCols(0) * s + j] / tss->getNumScenarios();
				dist = fabs(expval - floor(expval + 0.5));
				if (dist > maxdist) {
					maxdist = dist;
					branchingIndex = j;
					branchingValue = primsol[j];
				}
			}
		}
#endif
		findPhase++;
	}

	if (branchingIndex > -1) {

		branched = true;

		/** get branching index in first stage */
		if (branchingIndex < ncols_first_stage)
			branchingFirstStage = branchingIndex % tss->getNumCols(0);

		/** creating branching objects */
		branchingUp = new DspBranchObj();
		branchingDn = new DspBranchObj();
		for (int j = 0; j < master->ncols_orig_; ++j) {
			if (master->ctype_orig_[j] == 'C') continue;
			/** NOTE: branching on all the first-stage variables if SMIP */
			if (branchingIndex == j || (tss != NULL && branchingFirstStage == j % tss->getNumCols(0))) {
				DSPdebugMessage("Creating branch objects on column %d (value %e): [%e,%e] and [%e,%e]\n", 
					j, branchingValue, ceil(branchingValue), master->cubd_node_[j], master->clbd_node_[j], floor(branchingValue));
				branchingUp->push_back(j, ceil(branchingValue), master->cubd_node_[j]);
				branchingDn->push_back(j, master->clbd_node_[j], floor(branchingValue));
			} else if (master->clbd_node_[j] > master->clbd_orig_[j] || master->cubd_node_[j] < master->cubd_orig_[j]) {
				/** store any bound changes made in parent nodes */
				DSPdebugMessage("Adjusting bound change on column %d: [%e,%e]\n", j, master->clbd_node_[j], master->cubd_node_[j]);
				branchingUp->push_back(j, master->clbd_node_[j], master->cubd_node_[j]);
				branchingDn->push_back(j, master->clbd_node_[j], master->cubd_node_[j]);
			}
		}

		/** set best dual bounds */
		branchingUp->bestBound_ = master->getBestDualObjective();
		branchingDn->bestBound_ = master->getBestDualObjective();

		/** assign best dual solutions */
		branchingUp->dualsol_.assign(master->getBestDualSolution(), master->getBestDualSolution() + master->nrows_);
		branchingDn->dualsol_.assign(master->getBestDualSolution(), master->getBestDualSolution() + master->nrows_);

		/** set branching directions */
		branchingUp->direction_ = 1;
		branchingDn->direction_ = -1;

		/** add branching objects */
		branchingObjs.push_back(branchingUp);
		branchingObjs.push_back(branchingDn);
	} else {
		DSPdebugMessage("No branch object is found.\n");
	}

	END_TRY_CATCH_RTN(;,false)

	return branched;
}
