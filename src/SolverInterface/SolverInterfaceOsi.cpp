/*
 * SolverInterfaceOsi.cpp
 *
 *  Created on: Feb 9, 2015
 *      Author: kibaekkim
 */

#include <assert.h>
#include <Utility/DspMacros.h>
#include <Utility/DspMessage.h>
#include "SolverInterface/SolverInterfaceOsi.h"

#include "CoinWarmStartBasis.hpp"

SolverInterfaceOsi::SolverInterfaceOsi(DspParams * par) :
	SolverInterface(par),
	si_(NULL),
	ws_(NULL)
{
	BGN_TRY_CATCH
	DSP_RTN_CHECK_THROW(initialize());
	END_TRY_CATCH(;)
}

/** copy constructor */
SolverInterfaceOsi::SolverInterfaceOsi(SolverInterfaceOsi * si) :
	SolverInterface(si->par_)
{
	si_ = si->getOSI()->clone();
	ws_ = si->getWarmStart()->clone();
}

/** copy constructor */
SolverInterfaceOsi::SolverInterfaceOsi(DspParams * par, OsiSolverInterface * si) :
	SolverInterface(par)
{
	si_ = si->clone();
	ws_ = si->getWarmStart()->clone();
}

SolverInterfaceOsi::~SolverInterfaceOsi()
{
	BGN_TRY_CATCH
	DSP_RTN_CHECK_THROW(finalize());
	END_TRY_CATCH(;)
}

/** clone */
SolverInterface * SolverInterfaceOsi::clone()
{
	return new SolverInterfaceOsi(this);
}

/** initialize solver interface */
DSP_RTN_CODE SolverInterfaceOsi::initialize()
{
	return DSP_RTN_OK;
}

/** finalize solver interface */
DSP_RTN_CODE SolverInterfaceOsi::finalize()
{
	FREE_PTR(si_);
	FREE_PTR(ws_);
	return DSP_RTN_OK;
}

/** load problem */
void SolverInterfaceOsi::loadProblem(
		OsiSolverInterface * si,
		const char * probname)
{
	CoinPackedMatrix * mat = new CoinPackedMatrix(*si->getMatrixByCol());
	char * ctype = new char [si->getNumCols()];
	for (int j = 0; j < si->getNumCols(); ++j)
	{
		if (si->getColType()[j] == 0)
			ctype[j] = 'C';
		else if (si->getColType()[j] == 1)
			ctype[j] = 'B';
		else
			ctype[j] = 'I';
	}
	loadProblem(mat, si->getColLower(), si->getColUpper(), si->getObjCoefficients(),
			const_cast<const char*>(ctype), si->getRowLower(), si->getRowUpper(), probname);

	/** free memory */
	FREE_PTR(mat);
	FREE_ARRAY_PTR(ctype);
}

/** load problem */
void SolverInterfaceOsi::loadProblem(
		CoinPackedMatrix * mat,
		const double * collb,
		const double * colub,
		const double * obj,
		const char * ctype,
		const double * rowlb,
		const double * rowub,
		const char * probname)
{
	assert(si_);
	si_->loadProblem(*mat, collb, colub, obj, rowlb, rowub);
	for (int j = 0; j < mat->getNumCols(); ++j)
	{
		if (ctype[j] == 'C')
			si_->setContinuous(j);
		else
			si_->setInteger(j);
	}
}

/** add row */
void SolverInterfaceOsi::addRow(int size, const int * indices, const double * vals, double lb, double ub)
{
	assert(si_);
	si_->addRow(size, indices, vals, lb, ub);
}

/** add cuts and return the number of cuts applied */
int SolverInterfaceOsi::addCuts(OsiCuts cuts, double effectiveness)
{
	int numApplied  = 0;
#if 0
	OsiCuts oc;
	for (int i = 0; i < cuts.sizeCuts(); ++i)
	{
		OsiRowCut * rc = cuts.rowCutPtr(i);
		OsiColCut * cc = cuts.colCutPtr(i);
		if (rc)
		{
			if (rc->effectiveness() > effectiveness)
			{
				oc.insert(rc);
				numApplied++;
			}
		}
		else if (cc)
		{
			if (cc->effectiveness() > effectiveness)
			{
				oc.insert(cc);
				numApplied++;
			}
		}
	}
	si_->applyCuts(oc);
	oc.dumpCuts();
#else
	OsiSolverInterface::ApplyCutsReturnCode acRc = si_->applyCuts(cuts, effectiveness);
	numApplied = acRc.getNumApplied();
#endif

	if (ws_)
	{
		/** retrieve basis information */
		CoinWarmStartBasis * basis = dynamic_cast<CoinWarmStartBasis*>(ws_->clone());
		assert(basis != NULL);

		/** resize basis */
		basis->resize(si_->getNumRows(), si_->getNumCols());

		/** set basis info for new cuts added */
		for (int i = 0; i < numApplied; i++)
			basis->setArtifStatus(si_->getNumRows() - numApplied + i, CoinWarmStartBasis::basic);
		setWarmStart(basis);
	}

	return numApplied;
}

/** solve */
void SolverInterfaceOsi::solve()
{
	assert(si_);
	if (si_->getNumIntegers() == 0)
	{
		if (ws_ == NULL)
			si_->initialSolve();
		else
		{
			/** set warmstart basis */
			if (si_->setWarmStart(ws_) == false)
				printf("Failed to set warm-start information in SolverInterfaceOsi()\n");
			si_->resolve();
		}

		/** get warmstart */
		FREE_PTR(ws_);
		ws_ = si_->getWarmStart();
		assert(ws_);
	}
	else
		si_->branchAndBound();
}

/** solution status */
DSP_RTN_CODE SolverInterfaceOsi::getStatus()
{
	if (si_->isProvenOptimal())
		return DSP_STAT_OPTIMAL;
	if (si_->isProvenDualInfeasible())
		return DSP_STAT_DUAL_INFEASIBLE;
	if (si_->isProvenPrimalInfeasible())
		return DSP_STAT_PRIM_INFEASIBLE;
	if (si_->isPrimalObjectiveLimitReached())
		return DSP_STAT_LIM_PRIM_OBJ;
	if (si_->isDualObjectiveLimitReached())
		return DSP_STAT_LIM_DUAL_OBJ;
	if (si_->isIterationLimitReached())
		return DSP_STAT_LIM_ITERorTIME;
	if (si_->isAbandoned())
		return DSP_STAT_STOPPED_UNKNOWN;
	return DSP_STAT_UNKNOWN;
}

/** get global primal bound (upper bound in minimization) */
double SolverInterfaceOsi::getPrimalBound()
{
	return si_->getObjValue();
}

/** get dual bound (lower bound in minimization) */
double SolverInterfaceOsi::getDualBound()
{
	return si_->getObjValue();
}

/** get solution values */
const double * SolverInterfaceOsi::getSolution()
{
	return si_->getColSolution();
}

/** set objective coefficients */
void SolverInterfaceOsi::setObjCoef(double * obj)
{
	si_->setObjective(obj);
}

/** set iteration limit */
void SolverInterfaceOsi::setIterLimit(int limit)
{
	si_->setIntParam(OsiMaxNumIteration, limit);
}

/** set print out level */
void SolverInterfaceOsi::setPrintLevel(int level)
{
	assert(si_);
	si_->messageHandler()->setLogLevel(level);
	//si_->setLogLevel(level);
}

/** set warm start information */
void SolverInterfaceOsi::setWarmStart(CoinWarmStart * ws)
{
	FREE_PTR(ws_);
	ws_ = ws;
}

/** write MPS file */
void SolverInterfaceOsi::writeMps(const char * filename)
{
	assert(si_);
	si_->writeMps(filename);
}


