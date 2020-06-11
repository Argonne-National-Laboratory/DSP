/*
 * DwSolverSerial.cpp
 *
 *  Created on: Dec 5, 2016
 *      Author: kibaekkim
 */

#include "AlpsKnowledgeBrokerSerial.h"
#include "Model/TssModel.h"
#include "Solver/DantzigWolfe/DwSolverSerial.h"
#include "Solver/DantzigWolfe/DwMaster.h"
#include "Solver/DantzigWolfe/DwBundleDual.h"
#include "Solver/DantzigWolfe/DwBundleDualSmip.h"
#include "Solver/DantzigWolfe/DwWorker.h"

DwSolverSerial::DwSolverSerial(
		DecModel *   model,  /**< model pointer */
		DspParams *  par,    /**< parameters */
		DspMessage * message /**< message pointer */):
DecSolver(model, par, message),
master_(NULL),
worker_(NULL),
alps_(NULL) {

	/** just reset the files */
	std::ofstream ofs;
	if (par->getStrParam("VBC/FILE").size() > 0) {
		ofs.open(par->getStrParam("VBC/FILE").c_str());
		ofs.close();
	}
	if (par->getStrParam("DW/LOGFILE/OBJS").length() > 0) {
		ofs.open(par->getStrParam("DW/LOGFILE/OBJS").c_str());
		ofs.close();
	}

	/** time stamp */
	CoinWallclockTime();
}

DwSolverSerial::~DwSolverSerial() {
	FREE_PTR(master_);
	FREE_PTR(worker_);
	FREE_PTR(alps_);
}

DSP_RTN_CODE DwSolverSerial::init() {
	BGN_TRY_CATCH

	/** set parameter */
	par_->setIntPtrParamSize("ARR_PROC_IDX", model_->getNumSubproblems());
	for (int s = 0; s < model_->getNumSubproblems(); ++s)
		par_->setIntPtrParam("ARR_PROC_IDX", s, s);

	/** create worker */
	message_->print(1, "Initializing subproblems ... \n");
	worker_ = new DwWorker(model_, par_, message_);

	/** create master */
	message_->print(1, "Initializing master problem ... \n");
	//master_ = new DwMaster(worker_);
	if (model_->isStochastic())
		master_ = new DwBundleDualSmip(worker_);
	else
		master_ = new DwBundleDual(worker_);

	/** initialize master */
	DSP_RTN_CHECK_THROW(master_->init());

	/** create an Alps model */
	message_->print(1, "Initializing ALPS framework ... \n");
	if (model_->isStochastic())
		alps_ = new DwModelSmip(master_);
	else
		alps_ = new DwModel(master_);

	/** initialize the model */
	DSP_RTN_CHECK_THROW(alps_->init());

	/** parameter setting */
	DspParams* par = alps_->getSolver()->getParPtr();
	alps_->AlpsPar()->setEntry(AlpsParams::searchStrategy, par->getIntParam("ALPS/SEARCH_STRATEGY"));
	alps_->AlpsPar()->setEntry(AlpsParams::nodeLogInterval, par->getIntParam("ALPS/NODE_LOG_INTERVAL"));
	alps_->AlpsPar()->setEntry(AlpsParams::nodeLimit, par->getIntParam("ALPS/NODE_LIM"));
	alps_->AlpsPar()->setEntry(AlpsParams::timeLimit, par->getDblParam("ALPS/TIME_LIM"));
	alps_->AlpsPar()->setEntry(AlpsParams::clockType, AlpsClockTypeWallClock);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

DSP_RTN_CODE DwSolverSerial::solve() {
	BGN_TRY_CATCH

	/** solve */
	AlpsKnowledgeBrokerSerial alpsBroker(0, NULL, *alps_);
    alpsBroker.search(alps_);
	//alpsBroker.printBestSolution();
	
	AlpsExitStatus alpsstatus = alpsBroker.getSolStatus();
	switch(alpsstatus) {
	case AlpsExitStatusOptimal:
		status_ = DSP_STAT_OPTIMAL;
		break;
	case AlpsExitStatusTimeLimit:
	case AlpsExitStatusNodeLimit:
		status_ = DSP_STAT_LIM_ITERorTIME;
		break;
	case AlpsExitStatusSolLimit:
		status_ = DSP_STAT_STOPPED_SOLUTION;
		break;
	case AlpsExitStatusFeasible:
		status_ = DSP_STAT_FEASIBLE;
		break;
	case AlpsExitStatusInfeasible:
		status_ = DSP_STAT_PRIM_INFEASIBLE;
		break;
	case AlpsExitStatusFailed:
		status_ = DSP_STAT_ABORT;
		break;
	case AlpsExitStatusUnbounded:
		status_ = DSP_STAT_DUAL_INFEASIBLE;
		break;
	default:
		status_ = DSP_STAT_UNKNOWN;
		break;
	}

	DspNodeSolution* solution = NULL;
	if (alpsBroker.hasKnowledge(AlpsKnowledgeTypeSolution)) {
		solution = dynamic_cast<DspNodeSolution*>(alpsBroker.getBestKnowledge(AlpsKnowledgeTypeSolution).first);
		if (solution) {
			bestprimsol_ = solution->solution_;
			if (model_->isStochastic() && bestprimsol_.size() > 0) {
				TssModel* tss = dynamic_cast<TssModel*>(model_);
				if (bestprimsol_.size() > 0)
					bestprimsol_.erase(bestprimsol_.begin(), bestprimsol_.begin() + tss->getNumCols(0) * (tss->getNumScenarios() - 1));
			}
		}
	}
	bestprimobj_ = alpsBroker.getBestQuality();
	if (alpsstatus == AlpsExitStatusOptimal)
		bestdualobj_ = bestprimobj_;
	else
		bestdualobj_ = CoinMin(bestprimobj_, alps_->getBestDualObjective());

	message_->print(1, "Time spent in heuristics: %.2f seconds\n", alps_->heuristic_time_elapsed_);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

DSP_RTN_CODE DwSolverSerial::finalize() {
	BGN_TRY_CATCH

	/** finalize master */
	DSP_RTN_CHECK_THROW(master_->finalize());

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}
