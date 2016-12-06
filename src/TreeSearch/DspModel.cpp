/*
 * DspModel.cpp
 *
 *  Created on: Aug 31, 2016
 *      Author: kibaekkim
 */

/** Coin */
#include "AlpsEncoded.h"
#include "AlpsNodeDesc.h"
#include "AlpsKnowledgeBrokerSerial.h"
/** Dsp */
#include "Utility/DspMessage.h"
#include "TreeSearch/DspModel.h"
#include "TreeSearch/DspNodeDesc.h"
#include "TreeSearch/DspTreeNode.h"

DspModel::DspModel() :
		AlpsModel(),
		solver_(NULL) {
	/** nothing to do */
}

DspModel::DspModel(DecSolver* solver) :
		AlpsModel(),
		solver_(solver) {
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

	AlpsPar()->setEntry(AlpsParams::searchStrategy, 0);
	AlpsPar()->setEntry(AlpsParams::nodeLogInterval, 1);
	AlpsPar()->setEntry(AlpsParams::nodeLimit, 1);
	//AlpsPar()->setEntry(AlpsParams::timeLimit, 300);
	AlpsKnowledgeBrokerSerial alpsBroker(0, NULL, *this);
    alpsBroker.search(this);

//    DspNodeSolution* solution = dynamic_cast<DspNodeSolution*>(alpsBroker.getBestKnowledge(AlpsKnowledgeTypeSolution).first);
//    solution->print(std::cout);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}
