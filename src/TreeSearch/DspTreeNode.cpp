/*
 * DspTreeNode.cpp
 *
 *  Created on: Aug 30, 2016
 *      Author: kibaekkim
 */

// #define DSP_DEBUG

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
#include "Solver/DantzigWolfe/DwMaster.h"

DspTreeNode::DspTreeNode() : AlpsTreeNode() {}

DspTreeNode::~DspTreeNode() {
	for (auto obj = branchingObjs_.begin(); obj != branchingObjs_.end(); obj++) {
		FREE_PTR(*obj);
	}
	branchingObjs_.clear();
}

int DspTreeNode::process(bool isRoot, bool rampUp) {
	int status = AlpsReturnStatusOk;
	int ret; /**< solver return */

	/** retrieve objects */
	DspNodeDesc* desc = dynamic_cast<DspNodeDesc*>(desc_);
	DspModel* model = dynamic_cast<DspModel*>(desc_->getModel());
	DwMaster* solver = dynamic_cast<DwMaster*>(model->getSolver());
	DspParams* par = model->getParPtr();
	DspMessage* message = solver->getMessagePtr();

	/** bounds */
	double gUb = getKnowledgeBroker()->getIncumbentValue();
	double gLb = getKnowledgeBroker()->getBestNode()->getQuality();
	double parentLb = isRoot ? -ALPS_OBJ_MAX : getParent()->getQuality();
	double gap = (gUb - parentLb) / (fabs(gUb) + 1e-10);
	double relTol = par->getDblParam("DW/GAPTOL");
	DSPdebugMessage("Solving node %d, parentLb %e, gUb %e, gLb %e, gap %.2f\n", index_, parentLb, gUb, gLb, gap);

	/** fathom if the relative gap is small enough */
	if (gap < relTol) {
		//gap = (gUb - model->getBestDualObjective()) / (fabs(gUb) + 1e-10);
		model->setBestDualObjective(gLb);
		if (isCandidate()) {
			message->print(2, "The current node is pruned before processing (gap %.4f tol %.4f).\n", gap, relTol);
			setStatus(AlpsNodeStatusDiscarded);
			wirteLog("discarded", desc);
		} else {
			message->print(2, "The current node is fathomed (gap %.4f tol %.4f).\n", gap, relTol);
			setStatus(AlpsNodeStatusFathomed);
			wirteLog("fathomed", desc);
		}
		return status;
	}

	if (!isRoot) {
		/** set branching objects */
		model->setBranchingObjects(desc->getBranchingObject());

		if (par->getIntParam("DW/EVAL_UB") <= 0)
			par->setIntParam("DW/MAX_EVAL_UB", 0);
	}

	double alpsTimeRemain = par->getDblParam("ALPS/TIME_LIM") - getKnowledgeBroker()->timer().getWallClock();
	double dwTimeLim = CoinMin(par->getDblParam("DW/TIME_LIM"), alpsTimeRemain);
	model->setIterLimit(par->getIntParam("DW/ITER_LIM"));
	model->setTimeLimit(dwTimeLim);
	model->setBestPrimalObjective(std::min(gUb, ALPS_OBJ_MAX));

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
		gUb = model->getBestPrimalObjective();
		DSPdebugMessage("Found new upper bound %e\n", gUb);
		DspNodeSolution* nodesol = new DspNodeSolution(model->getBestPrimalSolution(), gUb);
		getKnowledgeBroker()->addKnowledge(AlpsKnowledgeTypeSolution, nodesol, gUb);
		wirteLog("heuristic", desc, gUb);
	}

	switch (model->getStatus()) {
	case DSP_STAT_OPTIMAL:
	case DSP_STAT_FEASIBLE:
	case DSP_STAT_LIM_DUAL_OBJ:
	case DSP_STAT_LIM_ITERorTIME: {
		// The current LB is obtained by the Dantzig-Wolfe problem (i.e., primal form).
		quality_ = model->getPrimalObjective();

		if (isRoot) {
			gLb = quality_;
			model->setBestDualObjective(gLb);
		}

//#define WRITE_PRIM_SOL
#ifdef WRITE_PRIM_SOL
		{
			char primsol_filename[128];
			sprintf(primsol_filename, "primsol%d.txt", index_);
			std::ofstream fp_primsol(primsol_filename);
			for (unsigned i = 0; i < solver->getLastSubprobSolutions().size(); ++i) {
				CoinPackedVector* sol = solver->getLastSubprobSolutions()[i];
				for (int j = 0; j < sol->getNumElements(); ++j)
					fp_primsol << sol->getIndices()[j] << "," << sol->getElements()[j] << std::endl;
				sol = NULL;
			}
			fp_primsol.close();
		}
#endif

		message->print(1, "[%f] curLb %.8e, bestUb %.8e, bestLb %.8e\n",
			getKnowledgeBroker()->timer().getWallClock(), quality_, gUb, gLb);

		// write objective values to log file
		if (par->getStrParam("DW/LOGFILE/OBJS").length() > 0) {
			log_dualobjs_.open(par->getStrParam("DW/LOGFILE/OBJS").c_str(), ios::app);
			if (isRoot) {
				for (unsigned i = 0; i < solver->log_time_.size(); ++i)
					log_dualobjs_ << solver->log_time_[i] << "," << solver->log_bestdual_bounds_[i] << "," << solver->log_bestprim_bounds_[i] << std::endl;
			} else {
				log_dualobjs_ << CoinWallclockTime() << "," << gLb << "," << gUb << std::endl;
			}
			log_dualobjs_.close();
		}

		/** fathom if LB is larger than UB. */
		if (quality_ >= gUb || quality_ >= ALPS_OBJ_MAX) {
			message->print(1, "The current node is fathomed.\n");
			setStatus(AlpsNodeStatusFathomed);
			wirteLog("fathomed", desc);
		} else {

			/** set global dual bound */
			// model->setBestDualObjective(gLb);

			/** Branching otherwise */
			message->print(3, "need to branch?");
			bool hasObjs = model->infeasibility() > model->feastol() ? model->chooseBranchingObjects(branchingObjs_) : false;
			message->print(3, " %s\n", hasObjs ? "yes" : "no");

			if (hasObjs) {
				/** set solution estimate; the lower the better */
				// solEstimate_ = branchingObjs_[0]->solEstimate_;
				//message->print(1, "solEstimate_ %e\n", solEstimate_);

				DSPdebugMessage("Branching on the current node.\n");
				setStatus(AlpsNodeStatusPregnant);
			} else {
				DSPdebugMessage("The current node has feasible solution.\n");
				if (quality_ < gUb) {
					DspNodeSolution* nodesol = new DspNodeSolution(model->getPrimalSolution(), quality_);
					getKnowledgeBroker()->addKnowledge(AlpsKnowledgeTypeSolution, nodesol, quality_);
				}
				/** no branching object is found; we are done! */
				setStatus(AlpsNodeStatusFathomed);
				wirteLog("integer", desc);
			}
		}
		break;
	}
	case DSP_STAT_PRIM_INFEASIBLE:
		//message->print(1, "The current node is infeasible.\n");
		setStatus(AlpsNodeStatusFathomed);
		wirteLog("infeasible", desc);
		break;
	default:
		message->print(1, "Unexpected solution status: %d.\n", model->getStatus());
		setStatus(AlpsNodeStatusDiscarded);
		wirteLog("fathomed", desc);
		//status = AlpsReturnStatusErr;
		break;
	}

	return status;
}

std::vector<CoinTriple<AlpsNodeDesc*, AlpsNodeStatus, double> > DspTreeNode::branch() {

	int ret;
	int solver_loglevel;

	/** retrieve objects */
	DspNodeDesc* desc = dynamic_cast<DspNodeDesc*>(desc_);
	DspModel* model = dynamic_cast<DspModel*>(desc->getModel());
	DwMaster* solver = dynamic_cast<DwMaster*>(model->getSolver());
	DspParams* par = model->getParPtr();
	DspMessage* message = solver->getMessagePtr();

	/** new nodes to be returned */
	std::vector<CoinTriple<AlpsNodeDesc*, AlpsNodeStatus, double> > newNodes;
	DspNodeDesc* node = NULL;

	/** bounds */
	double gUb = getKnowledgeBroker()->getIncumbentValue();
	double gap = (gUb - getQuality()) / (fabs(gUb) + 1e-10);
	double relTol = par->getDblParam("DW/GAPTOL");

	/** fathom if the relative gap is small enough */
	if (gap < relTol) {
		message->print(2, "The current node is fathomed (gap %.4f tol %.4f).\n", gap, relTol);
		setStatus(AlpsNodeStatusFathomed);
		wirteLog("fathomed", desc);
		return newNodes;
	}

	/** set status */
	setStatus(AlpsNodeStatusBranched);
	wirteLog("branched", desc, getQuality(), 1.0, 1);

	for (auto obj = branchingObjs_.begin(); obj != branchingObjs_.end(); obj++) {
		/** set branching object */
		model->setBranchingObjects(*obj);

		/** add branching-down node */
		node = new DspNodeDesc(model, (*obj)->direction_, *obj);

		/** TODO: Strong branching might not be a good idea for now,
			because there is no good implementation for warm-starting the node.
			Therefore, the strong branching is disabled by default. */

		/** Do strong down-branching */
		if (par->getBoolParam("DW/STRONG_BRANCH")) {
			message->print(2, "Strong branching ...\n");

			/** turn off display */
			solver_loglevel = solver->getLogLevel();
			solver->setLogLevel(0);
			/** set other parameters */
			bool run_heuristics = par->getBoolParam("DW/HEURISTICS");
			par->setBoolParam("DW/HEURISTICS", false);
			par->setBoolParam("DW/MASTER/REUSE_COLS", true);

			// solver->setIterLimit(10);
			ret = model->solve();

			/** restore solver display option */
			solver->setLogLevel(solver_loglevel);
			message->print(2, "ret %d status %d obj %e\n", ret, model->getStatus(), model->getDualObjective());
			par->setBoolParam("DW/HEURISTICS", run_heuristics);
			par->setBoolParam("DW/MASTER/REUSE_COLS", false);

			if (ret != DSP_RTN_OK) {
				continue;
			} else {
				gap = (gUb - getQuality()) / (fabs(gUb) + 1e-10);
				if (model->getStatus() == DSP_STAT_PRIM_INFEASIBLE ||
					model->getStatus() == DSP_STAT_LIM_DUAL_OBJ ||
					gap < relTol) {
					continue;
				} else {
					newNodes.push_back(CoinMakeTriple(
							static_cast<AlpsNodeDesc*>(node),
							AlpsNodeStatusCandidate,
							model->getDualObjective()));
					wirteLog("candidate", node, model->getDualObjective());
					//printf("Branching estimates objective value %e.\n", model->getDualObjective());
				}
			}
		} else {
			newNodes.push_back(CoinMakeTriple(
					static_cast<AlpsNodeDesc*>(node),
					AlpsNodeStatusCandidate,
					getQuality()));
			wirteLog("candidate", node, getQuality());
		}
		node = NULL;

	}

	return newNodes;
}

DspTreeNode* DspTreeNode::createNewTreeNode(AlpsNodeDesc*& desc) const {
	/** create a new node */
	DspTreeNode* node = new DspTreeNode();
	node->desc_ = desc;

	return node;
}

void DspTreeNode::wirteLog(const char* status, DspNodeDesc* desc, double lpbound, double infeas, int suminfeas) {

	DspModel* model = dynamic_cast<DspModel*>(desc->getModel());
	DspParams* par = model->getParPtr();

	if (par->getStrParam("VBC/FILE").size() > 0) {
		logstream_.open(par->getStrParam("VBC/FILE").c_str(), ios::app);

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

		logstream_.close();
	}
}
