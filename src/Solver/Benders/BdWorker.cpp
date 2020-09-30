/*
 * BdWorker.cpp
 *
 *  Created on: Feb 18, 2016
 *      Author: kibaekkim
 */

//#define DSP_DEBUG
#include "Solver/Benders/BdWorker.h"
#include "CoinHelperFunctions.hpp"

BdWorker::BdWorker(
		DecModel *   model,  /**< model pointer */
		DspParams *  par,    /**< parameters */
		DspMessage * message /**< message pointer */) :
model_(model),
par_(par),
message_(message) {

	/** parameters */
	parProcIdxSize_ = par_->getIntPtrParamSize("ARR_PROC_IDX");
	parProcIdx_     = par_->getIntPtrParam("ARR_PROC_IDX");

	/** create Benders subproblem solver */
	bdsub_ = new BdSub(par_);
	bdsub_->setSubIndices(parProcIdxSize_, parProcIdx_);
	bdsub_->loadProblem(model_);
}

BdWorker::BdWorker(const BdWorker& rhs) :
model_(rhs.model_),
par_(rhs.par_),
message_(rhs.message_),
parProcIdxSize_(rhs.parProcIdxSize_),
parProcIdx_(rhs.parProcIdx_) {
	bdsub_ = rhs.bdsub_->clone();
}

BdWorker::~BdWorker() {
	FREE_PTR(bdsub_);
}
