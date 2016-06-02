/*
 * SolverInterfaceSpx.cpp
 *
 *  Created on: Feb 9, 2015
 *      Author: kibaekkim
 */

#include <Utility/DspMacros.h>
#include "SolverInterface/SolverInterfaceSpx.h"
#include "OsiSpxSolverInterface.hpp"

SolverInterfaceSpx::SolverInterfaceSpx(DspParams * par) :
SolverInterfaceOsi(par)
{
	DSP_RTN_CHECK_THROW(initialize(), "initialize", "SolverInterfaceSpx");
}

/** copy constructor */
SolverInterfaceSpx::SolverInterfaceSpx(SolverInterfaceSpx * si) :
		SolverInterfaceOsi(si->par_)
{
	si_ = si->getOSI()->clone();
	ws_ = si->getWarmStart()->clone();
}

/** copy constructor */
SolverInterfaceSpx::SolverInterfaceSpx(DspParams * par, OsiSolverInterface * si) :
		SolverInterfaceOsi(par)
{
	OsiSpxSolverInterface * spx = dynamic_cast<OsiSpxSolverInterface*>(si);
	if (spx)
	{
		si_ = new OsiSpxSolverInterface;
		si_->loadProblem(*(si->getMatrixByCol()), si->getColLower(), si->getColUpper(),
				si->getObjCoefficients(), si->getRowLower(), si->getRowUpper());
		//si_ = si->clone();
	}
	ws_ = si->getWarmStart()->clone();
}

/** initialize solver interface */
DSP_RTN_CODE SolverInterfaceSpx::initialize()
{
	if (si_ == NULL)
	{
		si_ = new OsiSpxSolverInterface;
		si_->setHintParam(OsiDoDualInResolve);
		si_->setHintParam(OsiDoScale);
	}
	return DSP_RTN_OK;
}

