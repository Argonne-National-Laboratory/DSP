/*
 * SolverInterfaceCpx.cpp
 *
 *  Created on: Mar 14, 2016
 *      Author: kibaekkim
 */

/** Coin-OR */
#include "OsiCpxSolverInterface.hpp"

/** DSP */
#include "SolverInterface/SolverInterfaceCpx.h"

SolverInterfaceCpx::SolverInterfaceCpx(DspParams * par):
SolverInterfaceOsi(par), cpx_(NULL), useBarrier_(false)
{
	DSP_RTN_CHECK_THROW(initialize(), "initialize", "SolverInterfaceCpx");
}

/** copy constructor */
SolverInterfaceCpx::SolverInterfaceCpx(SolverInterfaceCpx * si) :
		SolverInterfaceOsi(si->par_), useBarrier_(si->useBarrier_)
{
	si_ = si->getOSI()->clone();
	ws_ = si->getWarmStart()->clone();
	cpx_ = dynamic_cast<OsiCpxSolverInterface*>(si_);
}

/** copy constructor */
SolverInterfaceCpx::SolverInterfaceCpx(DspParams * par, OsiSolverInterface * si) :
		SolverInterfaceOsi(par), useBarrier_(false)
{
	si_ = si->clone();
	ws_ = si->getWarmStart()->clone();
	cpx_ = dynamic_cast<OsiCpxSolverInterface*>(si_);
}

SolverInterfaceCpx::~SolverInterfaceCpx()
{
	cpx_ = NULL;
}

void SolverInterfaceCpx::solve()
{
	if (useBarrier_)
	{
		CPXbaropt(cpx_->getEnvironmentPtr(), cpx_->getLpPtr(OsiCpxSolverInterface::FREECACHED_RESULTS));
		int lpstat = CPXgetstat(cpx_->getEnvironmentPtr(), cpx_->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL));
		printf("CPLEX barrier status %d\n", lpstat);
	}
	else
		SolverInterfaceOsi::solve();
}

DSP_RTN_CODE SolverInterfaceCpx::initialize()
{
	if (si_ == NULL)
	{
		si_ = new OsiCpxSolverInterface;
		cpx_ = dynamic_cast<OsiCpxSolverInterface*>(si_);
		CPXsetintparam(cpx_->getEnvironmentPtr(), CPX_PARAM_THREADS, 1);
	}
	return DSP_RTN_OK;
}

void SolverInterfaceCpx::setNodeLimit(int limit)
{
	CPXsetintparam(cpx_->getEnvironmentPtr(), CPX_PARAM_NODELIM, limit);
}

void SolverInterfaceCpx::setIterLimit(int limit)
{
	CPXsetintparam(cpx_->getEnvironmentPtr(), CPX_PARAM_ITLIM, limit);
	CPXsetintparam(cpx_->getEnvironmentPtr(), CPX_PARAM_NETITLIM, limit);
	CPXsetintparam(cpx_->getEnvironmentPtr(), CPX_PARAM_STRONGITLIM, limit);
	CPXsetintparam(cpx_->getEnvironmentPtr(), CPX_PARAM_BARITLIM, limit);
}

void SolverInterfaceCpx::setTimeLimit(double sec)
{
	CPXsetdblparam(cpx_->getEnvironmentPtr(), CPX_PARAM_TILIM, CoinMin(sec,1.0e+75));
}

void SolverInterfaceCpx::setGapTol(double tol)
{
	CPXsetdblparam(cpx_->getEnvironmentPtr(), CPX_PARAM_EPGAP, CoinMin(tol,1.0));
}
