/*
 * DspTreeNode.cpp
 *
 *  Created on: Aug 30, 2016
 *      Author: kibaekkim
 */

#define DSP_DEBUG
/** Coin */
#include "CoinHelperFunctions.hpp"
#include "CoinUtility.hpp"
#include "AlpsKnowledgeBroker.h"
/** Dsp */
#include "Utility/DspMessage.h"
#include "Utility/DspRtnCodes.h"
#include "TreeSearch/DspTreeNode.h"
#include "TreeSearch/DspNodeDesc.h"
#include "TreeSearch/DspModel.h"
#include "TreeSearch/DspNodeSolution.h"

DspTreeNode::DspTreeNode() :
		AlpsTreeNode(),
		branchingUp_(NULL),
		branchingDn_(NULL) {
}

DspTreeNode::~DspTreeNode() {
	// TODO Auto-generated destructor stub
}

int DspTreeNode::process(bool isRoot, bool rampUp) {
	int status = AlpsReturnStatusOk;
	int ret; /**< solver return */
	double relTol = 0.0001;

	/** retrieve objects */
	DspNodeDesc* desc = dynamic_cast<DspNodeDesc*>(desc_);
	DspModel* model = dynamic_cast<DspModel*>(desc_->getModel());
	DecSolver* solver = model->getSolver();

	/** bounds */
	double curUb = getKnowledgeBroker()->getIncumbentValue();
	double gUb = curUb;
	double gLb = getKnowledgeBroker()->getBestNode()->getQuality();
	double gap = (gUb - gLb) / (fabs(gUb) + 1e-10);

	/** quality_ represents the best-known lower bound */
	if (isRoot)
		quality_ = -ALPS_OBJ_MAX;

#ifdef DSP_DEBUG
	if (depth_ > 100) {
		setStatus(AlpsNodeStatusFathomed);
		return status;
	}
#endif

	/** fathom if the relative gap is small enough */
	if (!isRoot) {
		if (gap < relTol) {
			setStatus(AlpsNodeStatusFathomed);
			return status;
		}
	}

	DSPdebugMessage("Solving node %d, gUb %e, gLb %e, gap %.2f\n", index_, gUb, gLb, gap);

	/** set branching object */
	solver->setBranchingObjects(desc->getBranchingObject());

	/** TODO: set bound information */
	//solver->setIterLimit(1);

	/** set heuristics */
	if (isRoot || gap > 0.1) {
		solver->setHeuristicRuns(true);
	} else {
		solver->setHeuristicRuns(false);
		DSPdebugMessage("Heuristics are turned off.\n");
	}

	while (1) {
		/** solve the bounding problem */
		ret = solver->solve();
		if (ret != DSP_RTN_OK) {
			setStatus(AlpsNodeStatusFathomed);
			return AlpsReturnStatusErr;
		}
		//DSPdebugMessage("Bounding solution status: %d\n", solver->getStatus());

		/** any heuristic solution */
		if (solver->getBestPrimalObjective() < gUb) {
			DSPdebugMessage("Found new upper bound %e\n", solver->getBestPrimalObjective());
			DspNodeSolution* nodesol = new DspNodeSolution(solver->getNumCols(),
					solver->getBestPrimalSolution(), solver->getBestPrimalObjective());
			getKnowledgeBroker()->addKnowledge(AlpsKnowledgeTypeSolution,
					nodesol, solver->getBestPrimalObjective());
		}

		if (solver->getStatus() != DSP_STAT_MW_RESOLVE)
			break;

		/** turn off heuristics for resolve */
		solver->setHeuristicRuns(false);
	}

	switch (solver->getStatus()) {
	case DSP_STAT_OPTIMAL:
	case DSP_STAT_FEASIBLE:{

		double curLb = solver->getDualObjective();
		double curUb = solver->getPrimalObjective();
		gUb = getKnowledgeBroker()->getIncumbentValue();

		/** update the best-known lower bound */
		if (curLb > quality_)
			quality_ = curLb;
		DSPdebugMessage("curLb %e, curUb %e, bestUb %e, bestLb %e\n", curLb, curUb, gUb, quality_);

		/** fathom if LB is larger than UB. */
		if (quality_ >= gUb)
			setStatus(AlpsNodeStatusFathomed);
		else {
			/** Branching otherwise */
			bool hasObjs = solver->chooseBranchingObjects(branchingUp_, branchingDn_);

			if (hasObjs) {
				DSPdebugMessage("Branching on the current node.\n");
				setStatus(AlpsNodeStatusPregnant);
			} else {
				DSPdebugMessage("The current node has feasible solution.\n");
				if (curUb < gUb) {
					DspNodeSolution* nodesol = new DspNodeSolution(
							solver->getNumCols(), solver->getPrimalSolution(), curUb);
					getKnowledgeBroker()->addKnowledge(AlpsKnowledgeTypeSolution, nodesol, curUb);
				}
				/** no branching object is found; we are done! */
				setStatus(AlpsNodeStatusFathomed);
			}
		}
		break;
	}
	case DSP_STAT_PRIM_INFEASIBLE:
		DSPdebugMessage("Fathom the current node.\n");
		setStatus(AlpsNodeStatusFathomed);
		break;
	default:
		assert(0);
		break;
	}

	return status;
}

std::vector<CoinTriple<AlpsNodeDesc*, AlpsNodeStatus, double> > DspTreeNode::branch() {

	/** retrieve objects */
	DspNodeDesc* desc = dynamic_cast<DspNodeDesc*>(desc_);
	DspModel* model = dynamic_cast<DspModel*>(desc->getModel());

	/** TODO: do strong branching */

	/** new nodes to be returned */
	std::vector<CoinTriple<AlpsNodeDesc*, AlpsNodeStatus, double> > newNodes;
	DspNodeDesc* node = NULL;

	/** add branching-up node */
	node = new DspNodeDesc(model, 1, branchingUp_);
	newNodes.push_back(CoinMakeTriple(
			static_cast<AlpsNodeDesc*>(node),
			AlpsNodeStatusCandidate,
			getQuality()));
	node = NULL;

	/** add branching-down node */
	node = new DspNodeDesc(model, -1, branchingDn_);
	newNodes.push_back(CoinMakeTriple(
			static_cast<AlpsNodeDesc*>(node),
			AlpsNodeStatusCandidate,
			getQuality()));
	node = NULL;

	/** set status */
	setStatus(AlpsNodeStatusBranched);

	return newNodes;
}

DspTreeNode* DspTreeNode::createNewTreeNode(AlpsNodeDesc*& desc) const {
	/** create a new node */
	DspTreeNode* node = new DspTreeNode();
	node->desc_ = desc;
	return node;
}
