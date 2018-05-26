/*
 * DwBranchNonant.cpp
 *
 *  Created on: May 2, 2018
 *      Author: Kibaek Kim
 */

// #define DSP_DEBUG
#include <DantzigWolfe/DwBranchNonant.h>

DwBranchNonant::DwBranchNonant(DwModel* model) : DwBranch(model) {
	DecSolver* solver = model_->getSolver();
	master_ = dynamic_cast<DwMaster*>(solver);
	if (solver->getModelPtr()->isStochastic())
		tss_ = dynamic_cast<TssModel*>(solver->getModelPtr());
}

bool DwBranchNonant::chooseBranchingObjects(
		std::vector<DspBranchObj*>& branchingObjs /**< [out] branching objects */) {

	if (tss_ == NULL) return false;

	bool branched = false;
	int branchingIndex = -1;
	double branchingValue, branchingDownValue, branchingUpValue;

	std::vector<double> refsol; // reference solution
	std::vector<double> devsol; // devations from the refsol
	std::vector<std::vector<double>> densesol; // dense solution vector

	DspBranchObj* branchingUp = NULL;
	DspBranchObj* branchingDn = NULL;

	BGN_TRY_CATCH

	/** retreive message object */
	DspMessage* message = model_->getSolver()->getMessagePtr();

	/** calculate reference solution */
	getRefSol(refsol);

	/** calculate variances from the reference solution */
	getDevSol(refsol, devsol);

#ifdef DSP_DEBUG
	printf("Reference Solution:\n");
	DspMessage::printArray(refsol.size(), &refsol[0]);
	printf("Deviations:\n");
	DspMessage::printArray(devsol.size(), &devsol[0]);
#endif

	double maxdev = 0.0;
	for (int j = 0; j < tss_->getNumCols(0); ++j) {
		if (devsol[j] > epsilon_) {
			maxdev = devsol[j];
			branchingIndex = j;
			branchingValue = refsol[j];
			if (tss_->getCtypeCore(0)[j] == 'C') {
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
			branchingIndex, branchingValue, master_->clbd_node_[branchingIndex], branchingDownValue, branchingUpValue, master_->cubd_node_[branchingIndex]);

		/** creating branching objects */
		branchingUp = new DspBranchObj();
		branchingDn = new DspBranchObj();
		for (int j = 0; j < tss_->getNumCols(0) * tss_->getNumScenarios(); ++j) {
			if (branchingIndex == j % tss_->getNumCols(0)) {
				DSPdebugMessage("Creating branch objects on column %d (value %e): [%e,%e] and [%e,%e]\n", 
					j, branchingValue, branchingUpValue, master_->cubd_node_[j], master_->clbd_node_[j], branchingDownValue);
				branchingUp->push_back(j, branchingUpValue, master_->cubd_node_[j]);
				branchingDn->push_back(j, master_->clbd_node_[j], branchingDownValue);
			} else if (master_->clbd_node_[j] > master_->clbd_orig_[j] || master_->cubd_node_[j] < master_->cubd_orig_[j]) {
				/** store any bound changes made in parent nodes */
				DSPdebugMessage("Adjusting bound change on column %d: [%e,%e]\n", j, master_->clbd_node_[j], master_->cubd_node_[j]);
				branchingUp->push_back(j, master_->clbd_node_[j], master_->cubd_node_[j]);
				branchingDn->push_back(j, master_->clbd_node_[j], master_->cubd_node_[j]);
			}
		}

		/** set best dual bounds */
		branchingUp->bestBound_ = master_->getBestDualObjective();
		branchingDn->bestBound_ = master_->getBestDualObjective();

		/** assign best dual solutions */
		branchingUp->dualsol_.assign(master_->getBestDualSolution(), master_->getBestDualSolution() + master_->nrows_);
		branchingDn->dualsol_.assign(master_->getBestDualSolution(), master_->getBestDualSolution() + master_->nrows_);

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

void DwBranchNonant::getRefSol(std::vector<double>& refsol) {
	refsol.resize(tss_->getNumCols(0), 0.0);
	for (unsigned s = 0; s < master_->getLastSubprobSolutions().size(); ++s) {
		const CoinPackedVector* sol = master_->getLastSubprobSolutions()[s];
		int sNumElements = sol->getNumElements();
		const int* sIndices = sol->getIndices();
		const double* sElements = sol->getElements();
		for (int j = 0; j < sNumElements; ++j)
			if (sIndices[j] < tss_->getNumCols(0) * tss_->getNumScenarios())
				refsol[sIndices[j] % tss_->getNumCols(0)] += sElements[j] * tss_->getProbability()[s];
	}
}

void DwBranchNonant::getDevSol(std::vector<double>& refsol, std::vector<double>& devsol) {
	devsol.resize(tss_->getNumCols(0), 0.0);
	for (unsigned s = 0; s < master_->getLastSubprobSolutions().size(); ++s) {
		const CoinPackedVector* sol = master_->getLastSubprobSolutions()[s];
		int sNumElements = sol->getNumElements();
		const int* sIndices = sol->getIndices();
		const double* sElements = sol->getElements();
		for (int j = 0; j < sNumElements; ++j)
			if (sIndices[j] < tss_->getNumCols(0) * tss_->getNumScenarios()) {
				int k = sIndices[j] % tss_->getNumCols(0);
				devsol[k] += pow(sElements[j] - refsol[k],2) * tss_->getProbability()[s];
			}
	}
}
