/*
 * DspModel.cpp
 *
 *  Created on: Aug 31, 2016
 *      Author: kibaekkim
 */

#include "AlpsEncoded.h"
#include "AlpsNodeDesc.h"
#include "Utility/DspMessage.h"
#include "TreeSearch/DspModel.h"
#include "TreeSearch/DspNodeDesc.h"
#include "TreeSearch/DspTreeNode.h"

DspModel::DspModel() :
		AlpsModel(),
		solver_(NULL), par_(NULL),
		status_(DSP_STAT_UNKNOWN),
		bestprimobj_(COIN_DBL_MAX),
		bestdualobj_(-COIN_DBL_MAX),
		primobj_(COIN_DBL_MAX),
		dualobj_(-COIN_DBL_MAX),
		infeasibility_(0.0) {
	/** nothing to do */
}

DspModel::DspModel(DecSolver* solver) :
		AlpsModel(),
		solver_(solver), par_(solver_->getParPtr()),
		status_(DSP_STAT_UNKNOWN),
		bestprimobj_(COIN_DBL_MAX),
		bestdualobj_(-COIN_DBL_MAX),
		primobj_(COIN_DBL_MAX),
		dualobj_(-COIN_DBL_MAX),
		infeasibility_(0.0) {
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
