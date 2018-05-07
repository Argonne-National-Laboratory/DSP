/*
 * DwBranchNonant.cpp
 *
 *  Created on: May 2, 2018
 *      Author: Kibaek Kim
 */

//#define DSP_DEBUG
#include <DantzigWolfe/DwBranchNonant.h>
#include <Model/TssModel.h>

bool DwBranchNonant::chooseBranchingObjects(
		std::vector<DspBranchObj*>& branchingObjs /**< [out] branching objects */) {

	DspMessage* message = model_->getSolver()->getMessagePtr();

	bool branched = false;
	int branchingIndex = -1;
	double branchingValue, branchingDownValue, branchingUpValue;

	DecSolver* solver = model_->getSolver();
	DwMaster* master = dynamic_cast<DwMaster*>(solver);
	int ncols_first_stage;
	std::vector<double> refsol; // reference solution
	std::vector<double> devsol; // devations from the refsol
	std::vector<std::vector<double>> densesol; // dense solution vector

	/** two-stage stochastic model */
	TssModel* tss = NULL;
	if (solver->getModelPtr()->isStochastic()) {
		tss = dynamic_cast<TssModel*>(solver->getModelPtr());
		ncols_first_stage = tss->getNumCols(0) * tss->getNumScenarios();
	} else {
		return branched;
	}

	if (master->getLastSubprobSolutions().size() != tss->getNumScenarios()) {
		message->print(0, "The number of subproblem solutions is not equal to the number of scenarios.\n");
		return branched;
	}

	DspBranchObj* branchingUp = NULL;
	DspBranchObj* branchingDn = NULL;

	BGN_TRY_CATCH

	/** calculate reference solution */
	refsol.resize(tss->getNumCols(0), 0.0);
	densesol.reserve(master->getLastSubprobSolutions().size());
	for (unsigned i = 0; i < master->getLastSubprobSolutions().size(); ++i) {
		CoinPackedVector* sol = master->getLastSubprobSolutions()[i];
#ifdef DSP_DEBUG2
		printf("Subproblem solution (scenario %d):\n", i);
		DspMessage::printArray(sol);
#endif
		std::vector<double> dsol(tss->getNumCols(0), 0.0);
		for (int j = 0; j < sol->getNumElements(); ++j) {
			if (sol->getIndices()[j] < ncols_first_stage) {
				refsol[sol->getIndices()[j] % tss->getNumCols(0)] += sol->getElements()[j] * tss->getProbability()[i];
				dsol[sol->getIndices()[j] % tss->getNumCols(0)] += sol->getElements()[j];
			}
		}
		densesol.push_back(dsol);
		sol = NULL;
	}

	/** calculate variances from the reference solution */
	devsol.resize(tss->getNumCols(0), 0.0);
	for (unsigned i = 0; i < master->getLastSubprobSolutions().size(); ++i) {
		CoinPackedVector* sol = master->getLastSubprobSolutions()[i];
		for (int j = 0; j < tss->getNumCols(0); ++j)
			devsol[j] += pow(densesol[i][j] - refsol[j],2) * tss->getProbability()[i];
		sol = NULL;
	}

#ifdef DSP_DEBUG
	printf("Reference Solution:\n");
	DspMessage::printArray(refsol.size(), &refsol[0]);
	printf("Deviations:\n");
	DspMessage::printArray(devsol.size(), &devsol[0]);
#endif

	double maxdev = 0.0;
	for (int j = 0; j < tss->getNumCols(0); ++j) {
		if (devsol[j] > epsilon_) {
			maxdev = devsol[j];
			branchingIndex = j;
			branchingValue = refsol[j];
			if (tss->getCtypeCore(0)[j] == 'C') {
				branchingDownValue = refsol[j] - epsilon_;
				branchingUpValue = refsol[j] + epsilon_;
			} else {
				branchingDownValue = floor(refsol[j]);
				branchingUpValue = ceil(refsol[j]);
			}
		}
	}

	if (branchingIndex > -1) {

		branched = true;

		message->print(2, "Creating branch objects on column %d (value %e): [%e,%e] and [%e,%e]\n", 
			branchingIndex, branchingValue, master->clbd_node_[branchingIndex], branchingDownValue, branchingUpValue, master->cubd_node_[branchingIndex]);

		/** creating branching objects */
		branchingUp = new DspBranchObj();
		branchingDn = new DspBranchObj();
		for (int j = 0; j < ncols_first_stage; ++j) {
			if (branchingIndex == j % tss->getNumCols(0)) {
				DSPdebugMessage("Creating branch objects on column %d (value %e): [%e,%e] and [%e,%e]\n", 
					j, branchingValue, branchingUpValue, master->cubd_node_[j], master->clbd_node_[j], branchingDownValue);
				branchingUp->push_back(j, branchingUpValue, master->cubd_node_[j]);
				branchingDn->push_back(j, master->clbd_node_[j], branchingDownValue);
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
