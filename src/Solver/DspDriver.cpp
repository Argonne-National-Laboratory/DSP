/*
 * DspDriver.cpp
 *
 *  Created on: Feb 16, 2016
 *      Author: kibaekkim
 */

#include "CoinFinite.hpp"
#include "Solver/DspDriver.h"

DspDriver::DspDriver(DspParams * par, DecModel * model):
	par_(par), model_(model),
	status_(DSP_STAT_UNKNOWN),
	primobj_(COIN_DBL_MAX), dualobj_(-COIN_DBL_MAX),
	cputime_(0.0), walltime_(0.0), numIterations_(0), numNodes_(0)
{
	/** create message */
	message_ = new DspMessage(par_->getIntParam("LOG_LEVEL"));
}

DspDriver::~DspDriver()
{
	par_ = NULL;
	model_ = NULL;
	FREE_PTR(message_);
	FREE_ARRAY_PTR(primsol_);
	FREE_ARRAY_PTR(dualsol_);
}

