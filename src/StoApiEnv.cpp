/*
 * StoApiEnv.cpp
 *
 *  Created on: Nov 6, 2014
 *      Author: kibaekkim
 */

#include "StoApiEnv.h"

StoApiEnv::StoApiEnv() : solver_(NULL), decdata_(NULL), model_(NULL)
{
	par_ = new StoParam;
}

StoApiEnv::~StoApiEnv()
{
	FREE_PTR(model_);
	FREE_PTR(solver_);
	FREE_PTR(par_);
	FREE_PTR(decdata_);
}

