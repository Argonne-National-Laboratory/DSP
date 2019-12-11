/*
 * DwBranchGenDisj.cpp
 *
 *  Created on: June 14, 2018
 *      Author: Kibaek Kim
 */

//#define DSP_DEBUG
#include "Solver/DantzigWolfe/DwBranchGenDisj.h"

DwBranchGenDisj::DwBranchGenDisj(DwModel* model) : DwBranch(model) {
	srand(0);
	DecSolver* solver = model_->getSolver();
	master_ = dynamic_cast<DwMaster*>(solver);
	if (solver->getModelPtr()->isStochastic())
		tss_ = dynamic_cast<TssModel*>(solver->getModelPtr());
}

bool DwBranchGenDisj::chooseBranchingObjects(
		std::vector<DspBranchObj*>& branchingObjs /**< [out] branching objects */) {

	if (tss_ == NULL) return false;

	bool branched = false;
	DspBranchObj* branchingUp = NULL;
	DspBranchObj* branchingDn = NULL;
	CoinPackedVector* vec = NULL;
	double rhs;

	BGN_TRY_CATCH

	/** retreive message object */
	DspMessage* message = model_->getSolver()->getMessagePtr();

	/** This example simply takes the sum of x per scenario 
		and creates a branch as follows:
			\sum_j x_j <= (average sum of x) - epsilon
		or
			\sum_j x_j >= (average sum of x) + epsilon

		Note: this is not a valid method.
	*/

	getSepHyperplane(&vec, rhs);

	double sumvec = 0.0;
	for (int i = 0; i < vec->getNumElements(); ++i)
		sumvec += vec->getElements()[i] * model_->getPrimalSolution()[vec->getIndices()[i]];
	double deviation = fabs(rhs - sumvec);
	printf("rhs %e sumvec %e deviation %e\n", rhs, sumvec, deviation);

	if (deviation > 1.0e-6) {

		branched = true;

		message->print(2, "Creating branch disjunction at value %e\n", rhs);

		branchingUp = new DspBranchObj();
		branchingUp->push_back(new CoinPackedVector(*vec), rhs + epsilon_, COIN_DBL_MAX);

		branchingDn = new DspBranchObj();
		branchingDn->push_back(new CoinPackedVector(*vec), -COIN_DBL_MAX, rhs - epsilon_);

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
		message->print(2, "No branch object is found.\n");
	}

	FREE_PTR(vec);

	END_TRY_CATCH_RTN(;,false)

	return branched;
}

void DwBranchGenDisj::getSepHyperplane(CoinPackedVector** vec, double& rhs) {
	*vec = new CoinPackedVector;
	(*vec)->reserve(tss_->getNumCols(0) * tss_->getNumScenarios());
	for (int j = 0; j < tss_->getNumCols(0); ++j) {
		double coef = (rand() % 100) / 100.0 + 0.5;
		(*vec)->insert(j, coef);
	}

	rhs = 0.0;
	for (int s = 0; s < tss_->getNumScenarios(); ++s)
		for (int j = 0; j < tss_->getNumCols(0); ++j)
			rhs += model_->getPrimalSolution()[s * tss_->getNumCols(0) + j] * tss_->getProbability()[s];
}
