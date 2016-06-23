/*
 * SolverInterfaceClp.cpp
 *
 *  Created on: Dec 9, 2014
 *      Author: kibaekkim
 */

#include <Utility/DspMacros.h>
#include "SolverInterface/SolverInterfaceClp.h"
#include "OsiClpSolverInterface.hpp"

SolverInterfaceClp::SolverInterfaceClp(DspParams * par) :
	SolverInterfaceOsi(par)
{
	DSP_RTN_CHECK_THROW(initialize(), "initialize", "SolverInterfaceClp");
}

/** copy constructor */
SolverInterfaceClp::SolverInterfaceClp(SolverInterfaceClp * si) :
		SolverInterfaceOsi(si->par_)
{
	si_ = si->getOSI()->clone();
	ws_ = si->getWarmStart()->clone();
}

/** copy constructor */
SolverInterfaceClp::SolverInterfaceClp(DspParams * par, OsiSolverInterface * si) :
		SolverInterfaceOsi(par)
{
	si_ = si->clone();
	ws_ = si->getWarmStart()->clone();
}

/** initialize solver interface */
DSP_RTN_CODE SolverInterfaceClp::initialize()
{
	if (si_ == NULL)
	{
		si_ = new OsiClpSolverInterface;
		si_->setHintParam(OsiDoDualInResolve);
		si_->setHintParam(OsiDoScale);
	}
	return DSP_RTN_OK;
}

