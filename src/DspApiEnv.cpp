/*
 * StoApiEnv.cpp
 *
 *  Created on: Nov 6, 2014
 *      Author: kibaekkim
 */

#include "DspApiEnv.h"
#include "Utility/DspMacros.h"

DspApiEnv::DspApiEnv() :
solver_(NULL), model_(NULL) {
	par_ = new DspParams;
	message_ = new DspMessage(par_->getIntParam("LOG_LEVEL"));
}

DspApiEnv::~DspApiEnv() {
	FREE_PTR(model_);
	FREE_PTR(solver_);
	FREE_PTR(par_);
	FREE_PTR(message_);
}

