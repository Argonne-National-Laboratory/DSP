/*
 * DspDriver.cpp
 *
 *  Created on: Feb 16, 2016
 *      Author: kibaekkim
 */

#include "DspDriver.h"
#include "Solver/DantzigWolfe/DwSolverSerial.h"

DspDriver::DspDriver(DecModel * model, DspParams * par):
		model_(model),
		par_(par),
		solver_(NULL),
		alps_(NULL),
		status_(DSP_STAT_UNKNOWN),
		primsol_(NULL),
		dualsol_(NULL),
		primobj_(COIN_DBL_MAX),
		dualobj_(-COIN_DBL_MAX),
		cputime_(0.0),
		walltime_(0.0),
		numIterations_(0),
		numNodes_(0) {
	/** create message */
	message_ = new DspMessage(par_->getIntParam("LOG_LEVEL"));
}

DspDriver::~DspDriver() {
	model_ = NULL;
	par_ = NULL;
	FREE_PTR(solver_);
	FREE_PTR(alps_);
	FREE_PTR(message_);
	FREE_ARRAY_PTR(primsol_);
	FREE_ARRAY_PTR(dualsol_);
}

DSP_RTN_CODE DspDriver::init() {
	BGN_TRY_CATCH

	solver_ = new DwSolverSerial(model_, par_, message_);
	DSP_RTN_CHECK_RTN_CODE(solver_->init());

	/** create an Alps model */
	alps_ = new DspModel(solver_);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DspDriver::run() {
	BGN_TRY_CATCH
	DSP_RTN_CHECK_RTN_CODE(alps_->solve());
	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

DSP_RTN_CODE DspDriver::finalize() {
	BGN_TRY_CATCH
	DSP_RTN_CHECK_RTN_CODE(solver_->finalize());
	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}
