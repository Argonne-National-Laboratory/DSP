/*
 * StoApiEnv.cpp
 *
 *  Created on: Nov 6, 2014
 *      Author: kibaekkim
 */

#include "DspApiEnv.h"
#include "DspMacros.h"

DspApiEnv::DspApiEnv() :
solver_(NULL), model_(NULL) {
	par_ = new DspParams;
}

DspApiEnv::~DspApiEnv() {
	FREE_PTR(model_);
	FREE_PTR(solver_);
	FREE_PTR(par_);
}

