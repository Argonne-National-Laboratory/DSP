/*
 * DspModel.cpp
 *
 *  Created on: Aug 31, 2016
 *      Author: kibaekkim
 */

/** Coin */
#include "AlpsEncoded.h"
#include "AlpsNodeDesc.h"
/** Dsp */
#include "Utility/DspMessage.h"
#include "TreeSearch/DspModel.h"
#include "TreeSearch/DspNodeDesc.h"
#include "TreeSearch/DspTreeNode.h"

DspModel::DspModel() :
		AlpsModel(),
		solver_(NULL), par_(NULL) {
	/** nothing to do */
}

DspModel::DspModel(DecSolver* solver) :
		AlpsModel(),
		solver_(solver), par_(solver_->getParPtr()) {
	/** nothing to do */
}

DspModel::~DspModel() {
	/** nothing to do */
}

AlpsTreeNode* DspModel::createRoot() {
	/** create root node */
	AlpsTreeNode* node = new DspTreeNode();

	/** create and set node description */
	DspNodeDesc* desc = new DspNodeDesc(this);
	node->setDesc(desc);

	return node;
}

bool DspModel::fathomAllNodes() {
	DSPdebugMessage("Fathom all nodes?\n");
	return false;
}

DSP_RTN_CODE DspModel::solve() {
	BGN_TRY_CATCH

	solver_->solve();

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

bool DspModel::chooseBranchingObjects(
		DspBranch*& branchingUp, /**< [out] branching-up object */
		DspBranch*& branchingDn  /**< [out] branching-down object */) {
#if 1
	return solver_->chooseBranchingObjects(branchingUp, branchingDn);
#else
	int findPhase = 0;
	bool branched = false;
	double dist, maxdist = 1.0e-6;
	int branchingIndex = -1;
	double branchingValue;

	BGN_TRY_CATCH

	/** cleanup */
	FREE_PTR(branchingUp)
	FREE_PTR(branchingDn)

	findPhase = 0;
	while (findPhase < 2 && branchingIndex < 0) {
		/** most fractional value */
		for (int j = 0; j < ncols_orig_; ++j) {
			if (ctype_orig_[j] == 'C') continue;
			dist = fabs(primsol_[j] - floor(primsol_[j] + 0.5));
			if (dist > maxdist) {
				maxdist = dist;
				branchingIndex = j;
				branchingValue = primsol_[j];
			}
		}
		findPhase++;
	}

	if (branchingIndex > -1) {
		DSPdebugMessage("Creating branch objects on column %d (value %e).\n", branchingIndex, branchingValue);
		branched = true;

		/** creating branching objects */
		branchingUp = new DspBranch();
		branchingDn = new DspBranch();
		for (int j = 0; j < ncols_orig_; ++j) {
			if (ctype_orig_[j] == 'C') continue;
			if (branchingIndex == j) {
				branchingUp->push_back(j, ceil(branchingValue), cubd_node_[j]);
				branchingDn->push_back(j, clbd_node_[j], floor(branchingValue));
			} else if (clbd_node_[j] > clbd_orig_[j] || cubd_node_[j] < cubd_orig_[j]) {
				/** store any bound changes made in parent nodes */
				branchingUp->push_back(j, clbd_node_[j], cubd_node_[j]);
				branchingDn->push_back(j, clbd_node_[j], cubd_node_[j]);
			}
		}
		branchingUp->bestBound_ = bestdualobj_;
		branchingDn->bestBound_ = bestdualobj_;
	} else {
		DSPdebugMessage("No branch object is found.\n");
	}

	END_TRY_CATCH_RTN(;,false)

	return branched;
#endif
}
