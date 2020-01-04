/*
 * DdWorkerCG.cpp
 *
 *  Created on: May 23, 2016
 *      Author: kibaekkim
 */

#include "DdWorkerCG.h"

DdWorkerCG::DdWorkerCG(
		DecModel *   model,  /**< model pointer */
		DspParams *  par,    /**< parameter pointer */
		DspMessage * message /**< message pointer */):
DdWorker(model, par, message) {
	/** parameters */
	parProcIdxSize_ = par_->getIntPtrParamSize("ARR_PROC_IDX");
	parProcIdx_     = par_->getIntPtrParam("ARR_PROC_IDX");
	DSPdebugMessage("parProcIdxSize_ %d\n", parProcIdxSize_);
}

DdWorkerCG::DdWorkerCG(const DdWorkerCG& rhs) :
DdWorker(rhs),
parProcIdxSize_(rhs.parProcIdxSize_),
parProcIdx_(rhs.parProcIdx_) {}
