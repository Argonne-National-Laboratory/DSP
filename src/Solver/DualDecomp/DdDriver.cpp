/*
 * DdDriver.cpp
 *
 *  Created on: Feb 16, 2016
 *      Author: kibaekkim
 */

//#define DSP_DEBUG

#include "Solver/DualDecomp/DdDriver.h"

/** constructor */
DdDriver::DdDriver(
		DspParams * par,
		DecModel * model):
DspDriver(par, model),
mw_(NULL)
{
	/** Nothing */
}

DdDriver::~DdDriver()
{
	FREE_PTR(mw_);
}
