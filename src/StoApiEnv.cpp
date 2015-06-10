/*
 * StoApiEnv.cpp
 *
 *  Created on: Nov 6, 2014
 *      Author: kibaekkim
 */

#include "StoApiEnv.h"

StoApiEnv::StoApiEnv() : solver_(NULL)
{
	par_ = new StoParam;
	tss_ = new TssModel;
}

StoApiEnv::~StoApiEnv()
{
	FREE_PTR(tss_);
	FREE_PTR(solver_);
	FREE_PTR(par_);
}

