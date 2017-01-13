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
	double relTol = 0.0001;

	/** retrieve objects */
	DspNodeDesc* desc = dynamic_cast<DspNodeDesc*>(desc_);
	DspModel* model = dynamic_cast<DspModel*>(desc_->getModel());
	DecSolver* solver = model->getSolver();
	DspParams* par = solver->getParPtr();

	double alpsTimeRemain = par->getDblParam("ALPS/TIME_LIM") - getKnowledgeBroker()->timer().getWallClock();
	double dwTimeLim = CoinMin(par->getDblParam("DW/TIME_LIM"), alpsTimeRemain);

	/** bounds */
	double curUb = getKnowledgeBroker()->getIncumbentValue();
	double gUb = curUb;
	double gLb = getKnowledgeBroker()->getBestNode()->getQuality();
	double gap = (gUb - gLb) / (fabs(gUb) + 1e-10);
	//printf("Solving node %d, gUb %e, gLb %e, gap %.2f\n", index_, gUb, gLb, gap);

	if (isRoot) {
		/** quality_ represents the best-known lower bound */
		quality_ = -ALPS_OBJ_MAX;

		/** set heuristics */
		solver->setHeuristicRuns(par->getBoolParam("DW/HEURISTICS"));
	} else {
		/** fathom if the relative gap is small enough */
		if (gap < relTol) {
			setStatus(AlpsNodeStatusFathomed);
			wirteLog("fathomed", desc);
			return status;
		}

		/** set heuristics */
		if (gap > 0.01)
			solver->setHeuristicRuns(par->getBoolParam("DW/HEURISTICS"));
		else
			solver->setHeuristicRuns(false);
		solver->setBranchingObjects(desc->getBranchingObject());
	}

	if (solver->getHeuristicRuns()) {
		/** The very first run takes one iteration and explores primal solutions. */
		solver->setIterLimit(1);
		solver->setTimeLimit(dwTimeLim);

		/** solve the bounding problem */
		ret = solver->solve();
		if (ret != DSP_RTN_OK) {
			setStatus(AlpsNodeStatusDiscarded);
			wirteLog("fathomed", desc);
			return AlpsReturnStatusErr;
		}
		DSPdebugMessage("Bounding solution status: %d\n", solver->getStatus());

		/** any heuristic solution */
		if (solver->getBestPrimalObjective() < gUb) {
			DSPdebugMessage("Found new upper bound %e\n", solver->getBestPrimalObjective());
			DspNodeSolution* nodesol = new DspNodeSolution(solver->getNumCols(),
					solver->getBestPrimalSolution(), solver->getBestPrimalObjective());
			getKnowledgeBroker()->addKnowledge(AlpsKnowledgeTypeSolution,
					nodesol, solver->getBestPrimalObjective());
			wirteLog("heuristic", desc, solver->getBestPrimalObjective());
		}
	}

	alpsTimeRemain = par->getDblParam("ALPS/TIME_LIM") - getKnowledgeBroker()->timer().getWallClock();
	dwTimeLim = CoinMin(par->getDblParam("DW/TIME_LIM"), alpsTimeRemain);

	solver->setIterLimit(par->getIntParam("DW/ITER_LIM"));
	solver->setTimeLimit(dwTimeLim);

	/** solve the bounding problem */
	ret = solver->solve();
	if (ret != DSP_RTN_OK) {
		setStatus(AlpsNodeStatusDiscarded);
		wirteLog("fathomed", desc);
		return AlpsReturnStatusErr;
	}
	DSPdebugMessage("Bounding solution status: %d\n", solver->getStatus());

	/** any heuristic solution */
	if (solver->getBestPrimalObjective() < gUb) {
		DSPdebugMessage("Found new upper bound %e\n", solver->getBestPrimalObjective());
		DspNodeSolution* nodesol = new DspNodeSolution(solver->getNumCols(),
				solver->getBestPrimalSolution(), solver->getBestPrimalObjective());
		getKnowledgeBroker()->addKnowledge(AlpsKnowledgeTypeSolution,
				nodesol, solver->getBestPrimalObjective());
		wirteLog("heuristic", desc, solver->getBestPrimalObjective());
	}

	switch (solver->getStatus()) {
	case DSP_STAT_OPTIMAL:
	case DSP_STAT_FEASIBLE:
	case DSP_STAT_LIM_ITERorTIME:
	{

		quality_ = solver->getDualObjective();
		double curUb = solver->getPrimalObjective();

		DSPdebugMessage("[%f] curLb %e, curUb %e, bestUb %e, bestLb %e\n",
				getKnowledgeBroker()->timer().getWallClock(), quality_, curUb, gUb, gLb);

		if (curUb - quality_ < -1.0e-8) {
			setStatus(AlpsNodeStatusDiscarded);
			wirteLog("fathomed", desc);
			return AlpsReturnStatusErr;
		}

		/** fathom if LB is larger than UB. */
		if (quality_ >= gUb) {
			setStatus(AlpsNodeStatusFathomed);
			wirteLog("fathomed", desc);
		} else {
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
		DSPdebugMessage("Unexpected solution status: %d.\n", solver->getStatus());
		assert(0);
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
	DecSolver* solver = model->getSolver();

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
	else if (strcmp(status, "candidate") == 0 || strcmp(status, "heuristic") == 0)
		logstream_ << " " << lpbound;
	logstream_ << std::endl;
#endif
}
