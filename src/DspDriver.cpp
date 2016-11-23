/*
 * DspDriver.cpp
 *
 *  Created on: Feb 16, 2016
 *      Author: kibaekkim
 */

#include "DspDriver.h"
#include "Solver/DantzigWolfe/DwMaster.h"
#include "Solver/DantzigWolfe/DwMasterTr.h"

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

	/** TODO: Parameterize what to create */
	solver_ = new DwMasterTr(model_, par_, message_);
	DSP_RTN_CHECK_RTN_CODE(solver_->init());

	alps_ = new DspModel(solver_);

	return DSP_RTN_OK;
}

DSP_RTN_CODE DspDriver::run() {
	DSP_RTN_CHECK_RTN_CODE(alps_->solve());
	//solver_->solve();
	return DSP_RTN_OK;
}

DSP_RTN_CODE DspDriver::finalize() {
	solver_->finalize();
	return DSP_RTN_OK;
}
