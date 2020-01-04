/*
 * DdDriver.cpp
 *
 *  Created on: Feb 16, 2016
 *      Author: kibaekkim
 */

//#define DSP_DEBUG

#include "Solver/DualDecomp/DdDriver.h"

DdDriver::DdDriver(
		DecModel *   model,  /**< model pointer */
		DspParams *  par,    /**< parameters */
		DspMessage * message /**< message pointer */):
DecSolver(model, par, message),
mw_(NULL) {}

DdDriver::DdDriver(const DdDriver& rhs) :
DecSolver(rhs) {
	mw_ = rhs.mw_->clone();
}

DdDriver::~DdDriver() {
	FREE_PTR(mw_);
}
