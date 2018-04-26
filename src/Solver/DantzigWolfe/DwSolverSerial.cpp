/*
 * DwSolverSerial.cpp
 *
 *  Created on: Dec 5, 2016
 *      Author: kibaekkim
 */

#include "AlpsKnowledgeBrokerSerial.h"
#include <Model/TssModel.h>
#include <DantzigWolfe/DwSolverSerial.h>
#include <DantzigWolfe/DwMaster.h>
#include <DantzigWolfe/DwBundleDual.h>
#include <DantzigWolfe/DwWorker.h>

DwSolverSerial::DwSolverSerial(
		DecModel *   model,  /**< model pointer */
		DspParams *  par,    /**< parameters */
		DspMessage * message /**< message pointer */):
DecSolver(model, par, message),
master_(NULL),
worker_(NULL),
alps_(NULL) {}

DwSolverSerial::~DwSolverSerial() {
	FREE_PTR(master_);
	FREE_PTR(worker_);
	FREE_PTR(alps_);
}

DSP_RTN_CODE DwSolverSerial::init() {
	BGN_TRY_CATCH

	/** create worker */
	message_->print(1, "Initializing subproblems ... \n");
	worker_ = new DwWorker(model_, par_, message_);

	/** create master */
	message_->print(1, "Initializing master problem ... \n");
	//master_ = new DwMaster(worker_);
	master_ = new DwBundleDual(worker_);

	/** initialize master */
	DSP_RTN_CHECK_THROW(master_->init());

	/** create an Alps model */
	message_->print(1, "Initializing ALPS framework ... \n");
	alps_ = new DwModel(master_);

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
	// alpsBroker.printBestSolution();

	DspNodeSolution* solution = NULL;
	if (alpsBroker.hasKnowledge(AlpsKnowledgeTypeSolution)) {
		solution = dynamic_cast<DspNodeSolution*>(alpsBroker.getBestKnowledge(AlpsKnowledgeTypeSolution).first);
		bestprimsol_ = solution->solution_;
		if (model_->isStochastic()) {
			TssModel* tss = dynamic_cast<TssModel*>(model_);
			bestprimsol_.erase(bestprimsol_.begin(), bestprimsol_.begin() + tss->getNumCols(0) * (tss->getNumScenarios() - 1));
		}
	}
	bestprimobj_ = alpsBroker.getBestQuality();
	if (alpsBroker.getSolStatus() == AlpsExitStatusOptimal)
		bestdualobj_ = bestprimobj_;
	else
		bestdualobj_ = alpsBroker.getBestKnowledge(AlpsKnowledgeTypeNode).second;

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
