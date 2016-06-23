/*
 * DdWorkerCG.cpp
 *
 *  Created on: May 23, 2016
 *      Author: kibaekkim
 */

#include "DdWorkerCG.h"

DdWorkerCG::DdWorkerCG(DspParams * par, DecModel * model, DspMessage * message):
	DdWorker(par, model, message) {
	/** parameters */
	parProcIdxSize_ = par_->getIntPtrParamSize("ARR_PROC_IDX");
	parProcIdx_     = par_->getIntPtrParam("ARR_PROC_IDX");
	DSPdebugMessage("parProcIdxSize_ %d\n", parProcIdxSize_);
}

DdWorkerCG::~DdWorkerCG() {
	// TODO Auto-generated destructor stub
}
