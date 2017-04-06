/*
 * DspTreeNode.cpp
 *
 *  Created on: Aug 30, 2016
 *      Author: kibaekkim
 */

//#define DSP_DEBUG
//#define WRITELOG

/** Coin */
#include "CoinHelperFunctions.hpp"
#include "CoinUtility.hpp"
#include "AlpsKnowledgeBroker.h"
/** Dsp */
#include "Utility/DspMessage.h"
#include "Utility/DspRtnCodes.h"
#include "TreeSearch/DspTreeNode.h"
#include "TreeSearch/DspModel.h"
#include "TreeSearch/DspNodeSolution.h"

DspTreeNode::DspTreeNode() :
		AlpsTreeNode(),
		branchingUp_(NULL),
		branchingDn_(NULL) {
#ifdef WRITELOG
	//const char* logname = getKnowledgeBroker()->getModel()->AlpsPar()->entry(AlpsParams::logFile).c_str();
	DSPdebugMessage("Writing log file to DspAlps.vbc.\n");
	logstream_.open("DspAlps.vbc", ios::app);
#endif
}

DspTreeNode::~DspTreeNode() {
	FREE_PTR(branchingUp_);
	FREE_PTR(branchingDn_);
#ifdef WRITELOG
	logstream_.close();
#endif
}

int DspTreeNode::process(bool isRoot, bool rampUp) {
	int status = AlpsReturnStatusOk;
	int ret; /**< solver return */

	/** retrieve objects */
	DspNodeDesc* desc = dynamic_cast<DspNodeDesc*>(desc_);
	DspModel* model = dynamic_cast<DspModel*>(desc_->getModel());
	DspParams* par = model->getParPtr();
	double relTol = 0.0001;//par->getDblParam("DW/GAPTOL");

	/** bounds */
	double curUb = getKnowledgeBroker()->getIncumbentValue();
	double gUb = curUb;
	double gLb = getKnowledgeBroker()->getBestNode()->getQuality();
	double gap = (gUb - gLb) / (fabs(gUb) + 1e-10);
	//printf("Solving node %d, gUb %e, gLb %e, gap %.2f\n", index_, gUb, gLb, gap);

	if (isRoot) {
		/** quality_ represents the best-known lower bound */
		quality_ = -ALPS_OBJ_MAX;
	} else {
		/** fathom if the relative gap is small enough */
		if (gap < relTol) {
			setStatus(AlpsNodeStatusFathomed);
			wirteLog("fathomed", desc);
			return status;
		}
		/** set branching objects */
		model->setBranchingObjects(desc->getBranchingObject());
	}

	double alpsTimeRemain = par->getDblParam("ALPS/TIME_LIM") - getKnowledgeBroker()->timer().getWallClock();
	double dwTimeLim = CoinMin(par->getDblParam("DW/TIME_LIM"), alpsTimeRemain);
	model->setIterLimit(par->getIntParam("DW/ITER_LIM"));
	model->setTimeLimit(dwTimeLim);

	/** solve the bounding problem */
	ret = model->solve();
	if (ret != DSP_RTN_OK) {
		setStatus(AlpsNodeStatusDiscarded);
		wirteLog("fathomed", desc);
		return AlpsReturnStatusErr;
	}
	DSPdebugMessage("Bounding solution status: %d\n", model->getStatus());

	/** any heuristic solution */
	if (model->getBestPrimalObjective() < gUb) {
		DSPdebugMessage("Found new upper bound %e\n", model->getBestPrimalObjective());
		DspNodeSolution* nodesol = new DspNodeSolution(model->getBestPrimalSolution(), model->getBestPrimalObjective());
		getKnowledgeBroker()->addKnowledge(AlpsKnowledgeTypeSolution, nodesol, model->getBestPrimalObjective());
		wirteLog("heuristic", desc, model->getBestPrimalObjective());
	}

	switch (model->getStatus()) {
	case DSP_STAT_OPTIMAL:
	case DSP_STAT_FEASIBLE:
	case DSP_STAT_LIM_ITERorTIME: {
		quality_ = model->getDualObjective();
		double curUb = model->getPrimalObjective();

		printf("[%f] curLb %.8e, curUb %.8e, bestUb %.8e, bestLb %.8e\n",
			getKnowledgeBroker()->timer().getWallClock(), quality_, curUb, gUb, gLb);

		/** fathom if LB is larger than UB. */
		if (quality_ >= gUb || curUb >= 1.0e+20) {
			setStatus(AlpsNodeStatusFathomed);
			wirteLog("fathomed", desc);
		} else {
			/** Branching otherwise */
			bool hasObjs = model->chooseBranchingObjects(branchingUp_, branchingDn_);

			if (hasObjs) {
				DSPdebugMessage("Branching on the current node.\n");
				setStatus(AlpsNodeStatusPregnant);
			} else {
				DSPdebugMessage("The current node has feasible solution.\n");
				if (curUb < gUb) {
					DspNodeSolution* nodesol = new DspNodeSolution(model->getPrimalSolution(), curUb);
					getKnowledgeBroker()->addKnowledge(AlpsKnowledgeTypeSolution, nodesol, curUb);
				}
				/** no branching object is found; we are done! */
				setStatus(AlpsNodeStatusFathomed);
				wirteLog("integer", desc);
			}
		}
		break;
	}
	case DSP_STAT_PRIM_INFEASIBLE:
		DSPdebugMessage("Fathom the current node.\n");
		setStatus(AlpsNodeStatusFathomed);
		wirteLog("infeasible", desc);
		break;
	default:
		DSPdebugMessage("Unexpected solution status: %d.\n", model->getStatus());
		setStatus(AlpsNodeStatusDiscarded);
		wirteLog("fathomed", desc);
		//status = AlpsReturnStatusErr;
		break;
	}

	return status;
}

std::vector<CoinTriple<AlpsNodeDesc*, AlpsNodeStatus, double> > DspTreeNode::branch() {

	int ret, status;
	int solver_loglevel;

	/** retrieve objects */
	DspNodeDesc* desc = dynamic_cast<DspNodeDesc*>(desc_);
	DspModel* model = dynamic_cast<DspModel*>(desc->getModel());

	/** new nodes to be returned */
	std::vector<CoinTriple<AlpsNodeDesc*, AlpsNodeStatus, double> > newNodes;
	DspNodeDesc* node = NULL;

	/** set status */
	setStatus(AlpsNodeStatusBranched);
	wirteLog("branched", desc, getQuality(), 1.0, 1);
//#define STRONG_BRANCH
#ifdef STRONG_BRANCH
	/** turn off display */
	solver_loglevel = solver->getLogLevel();
	solver->setLogLevel(0);
	/** set other parameters */
	solver->setHeuristicRuns(false);
	solver->setIterLimit(10);

	/** Do strong down-branching */
	solver->setBranchingObjects(branchingDn_);
	ret = solver->solve();

	/** add branching-down node */
	node = new DspNodeDesc(model, -1, branchingDn_);
	if (ret != DSP_RTN_OK) {
		newNodes.push_back(CoinMakeTriple(
				static_cast<AlpsNodeDesc*>(node),
				AlpsNodeStatusDiscarded,
				ALPS_OBJ_MAX));
	} else {
		if (solver->getStatus() == DSP_STAT_PRIM_INFEASIBLE) {
			newNodes.push_back(CoinMakeTriple(
					static_cast<AlpsNodeDesc*>(node),
					AlpsNodeStatusFathomed,
					ALPS_OBJ_MAX));
			wirteLog("infeasible", node);
			DSPdebugMessage("Strong branching fathomed the child.\n");
		} else {
			newNodes.push_back(CoinMakeTriple(
					static_cast<AlpsNodeDesc*>(node),
					AlpsNodeStatusCandidate,
					solver->getPrimalObjective()));
			wirteLog("candidate", node, solver->getPrimalObjective());
			DSPdebugMessage("Strong branching estimates objective value %e.\n", solver->getPrimalObjective());
		}
	}
	node = NULL;

	/** Do strong UP-branching */
	solver->setBranchingObjects(branchingUp_);
	ret = solver->solve();

	/** add branching-UP node */
	node = new DspNodeDesc(model, 1, branchingUp_);
	if (ret != DSP_RTN_OK) {
		newNodes.push_back(CoinMakeTriple(
				static_cast<AlpsNodeDesc*>(node),
				AlpsNodeStatusDiscarded,
				ALPS_OBJ_MAX));
	} else {
		if (solver->getStatus() == DSP_STAT_PRIM_INFEASIBLE) {
			newNodes.push_back(CoinMakeTriple(
					static_cast<AlpsNodeDesc*>(node),
					AlpsNodeStatusFathomed,
					ALPS_OBJ_MAX));
			wirteLog("infeasible", node);
			DSPdebugMessage("Strong branching fathomed the child.\n");
		} else {
			newNodes.push_back(CoinMakeTriple(
					static_cast<AlpsNodeDesc*>(node),
					AlpsNodeStatusCandidate,
					solver->getPrimalObjective()));
			wirteLog("candidate", node, solver->getPrimalObjective());
			DSPdebugMessage("Strong branching estimates objective value %e.\n", solver->getPrimalObjective());
		}
	}
	node = NULL;

	/** restore solver display option */
	solver->setLogLevel(solver_loglevel);
#else
	/** add branching-down node */
	node = new DspNodeDesc(model, -1, branchingDn_);
	newNodes.push_back(CoinMakeTriple(
			static_cast<AlpsNodeDesc*>(node),
			AlpsNodeStatusCandidate, getQuality()));
	wirteLog("candidate", node, getQuality());
	node = NULL;

	/** add branching-UP node */
	node = new DspNodeDesc(model, 1, branchingUp_);
	newNodes.push_back(CoinMakeTriple(
			static_cast<AlpsNodeDesc*>(node),
			AlpsNodeStatusCandidate, getQuality()));
	wirteLog("candidate", node, getQuality());
	node = NULL;
#endif

	return newNodes;
}

DspTreeNode* DspTreeNode::createNewTreeNode(AlpsNodeDesc*& desc) const {
	/** create a new node */
	DspTreeNode* node = new DspTreeNode();
	node->desc_ = desc;

	return node;
}

void DspTreeNode::wirteLog(const char* status, DspNodeDesc* desc, double lpbound, double infeas, int suminfeas) {
#ifdef WRITELOG
	logstream_ << getKnowledgeBroker()->timer().getWallClock() << " " << status << " ";
	if (strcmp(status, "candidate") == 0) {
		int nextindex = getKnowledgeBroker()->getNextNodeIndex();
		if (desc->branchdir() > 0)
			logstream_ << nextindex << " " << index_ << " ";
		else
			logstream_ << nextindex + 1 << " " << index_ << " ";
	} else if (index_ == 0)
		logstream_ << index_ << " " << index_ << " ";
	else
		logstream_ << index_ << " " << parentIndex_ << " ";
	if (index_ == 0)
		logstream_ << "M";
	else if (desc->branchdir() > 0)
		logstream_ << "L";
	else
		logstream_ << "R";
	if (strcmp(status, "branched") == 0)
		logstream_ << " " << lpbound << " " << infeas << " " << suminfeas;
	else if (strcmp(status, "candidate") == 0 || strcmp(status, "heuristic") == 0 || strcmp(status, "integer") == 0)
		logstream_ << " " << lpbound;
	logstream_ << std::endl;
#endif
}
