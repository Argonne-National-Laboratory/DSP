/*
 * DeDriver.cpp
 *
 *  Created on: Feb 17, 2016
 *      Author: kibaekkim
 */

/** Coin */
#include "OsiCbcSolverInterface.hpp"
#include "OsiCpxSolverInterface.hpp"
/** Dsp */
#include "Model/TssModel.h"
#include "Solver/Deterministic/DeDriver.h"

DeDriver::DeDriver(DspParams * par, DecModel * model):
	DspDriver(par,model), si_(NULL) {}

DeDriver::~DeDriver()
{
	FREE_PTR(si_);
}

/** initilize */
DSP_RTN_CODE DeDriver::init()
{
	BGN_TRY_CATCH

	primsol_ = new double [model_->getFullModelNumCols()];

	END_TRY_CATCH(;)

	return DSP_RTN_OK;
}

/** run */
DSP_RTN_CODE DeDriver::run()
{
#define FREE_MEMORY       \
	FREE_PTR(mat)         \
	FREE_ARRAY_PTR(clbd)  \
	FREE_ARRAY_PTR(cubd)  \
	FREE_ARRAY_PTR(ctype) \
	FREE_ARRAY_PTR(obj)   \
	FREE_ARRAY_PTR(rlbd)  \
	FREE_ARRAY_PTR(rubd)

	assert(model_);

	/** model info */
	CoinPackedMatrix * mat = NULL;
	double * clbd   = NULL;
	double * cubd   = NULL;
	double * obj    = NULL;
	char *   ctype  = NULL;
	double * rlbd   = NULL;
	double * rubd   = NULL;

	BGN_TRY_CATCH

	double stime;

	/** get DE model */
	DSP_RTN_CHECK_THROW(model_->getFullModel(mat, clbd, cubd, ctype, obj, rlbd, rubd));

	int nIntegers = model_->getNumIntegers();

	if (model_->isStochastic())
	{
		TssModel * tssModel;
		try
		{
			tssModel = dynamic_cast<TssModel *>(model_);
		}
		catch (const std::bad_cast& e)
		{
			printf("Model claims to be stochastic when it is not");
			return DSP_RTN_ERR;
		}

		/** relax integrality? */
		if (par_->getBoolPtrParam("RELAX_INTEGRALITY")[0])
		{
			for (int j = 0; j < tssModel->getNumCols(0); ++j)
			{
				if (ctype[j] != 'C')
					nIntegers--;
				ctype[j] = 'C';
			}
		}
		if (par_->getBoolPtrParam("RELAX_INTEGRALITY")[1])
		{
			for (int j = 0; j < tssModel->getNumCols(1); ++j)
			{
				if (ctype[tssModel->getNumCols(0) + j] != 'C')
					nIntegers--;
			}
			CoinFillN(ctype + tssModel->getNumCols(0), tssModel->getNumScenarios() * tssModel->getNumCols(1), 'C');
		}
	}
	else
	{
		if (par_->getBoolPtrParam("RELAX_INTEGRALITY")[0] ||
				par_->getBoolPtrParam("RELAX_INTEGRALITY")[1])
		{
			for (int j = 0; j < mat->getNumCols(); j++)
			{
				if (ctype[j] != 'C')
					nIntegers--;
				ctype[j] = 'C';
			}
		}
	}

	/** create solver */
	si_ = new OsiCbcSolverInterface();
	/** print level */
	si_->messageHandler()->logLevel(par_->getIntParam("LOG_LEVEL"));
	/** load problem */
	si_->loadProblem(*mat, clbd, cubd, obj, ctype, rlbd, rubd);

	double time_limit = CoinMin(
			par_->getDblParam("DE/WALL_LIM"),
			par_->getDblParam("SCIP/TIME_LIM"));
	/** TODO: lines specific to CBC */
	OsiCbcSolverInterface* cbc = dynamic_cast<OsiCbcSolverInterface*>(si_);
	if (cbc) {
		/** time limit */
		cbc->setMaximumSeconds(time_limit);
		/** set node limit */
		cbc->setMaximumNodes(par_->getIntParam("NODE_LIM"));
	}

	/** tic */
	cputime_  = CoinCpuTime();
	walltime_ = CoinGetTimeOfDay();

	/** solve */
	si_->initialSolve();
	if (nIntegers > 0)
		si_->branchAndBound();

	/** toc */
	cputime_  = CoinCpuTime() - cputime_;
	walltime_ = CoinGetTimeOfDay() - walltime_;

	/** solution status */
	if (si_->isDualObjectiveLimitReached())
		status_ = DSP_STAT_LIM_ITERorTIME;
	else if (si_->isIterationLimitReached())
		status_ = DSP_STAT_LIM_ITERorTIME;
	else if (si_->isPrimalObjectiveLimitReached())
		status_ = DSP_STAT_LIM_PRIM_OBJ;
	else if (si_->isDualObjectiveLimitReached())
		status_ = DSP_STAT_LIM_DUAL_OBJ;
	else if (si_->isProvenPrimalInfeasible())
		status_ = DSP_STAT_PRIM_INFEASIBLE;
	else if (si_->isProvenDualInfeasible())
		status_ = DSP_STAT_DUAL_INFEASIBLE;
	else if (si_->isProvenOptimal())
		status_ = DSP_STAT_OPTIMAL;

	/** get solutions */
	if (status_ != DSP_STAT_PRIM_INFEASIBLE &&
			status_ != DSP_STAT_DUAL_INFEASIBLE) {
		/** objective bounds */
		primobj_ = si_->getObjValue();
		/** TODO: how to get dual bound? */
		dualobj_ = si_->getObjValue();
		/** solution */
		CoinCopyN(si_->getColSolution(), si_->getNumCols(), primsol_);
		/** statistics */
		numIterations_ = si_->getIterationCount();
		/** TODO: Cbc specific */
		if (cbc)
			numNodes_ = cbc->getNodeCount();
	}

	/** save memory */
	FREE_MEMORY

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	return DSP_RTN_OK;

#undef FREE_MEMORY
}

DSP_RTN_CODE DeDriver::finalize() {
	return DSP_RTN_OK;
}

void DeDriver::writeExtMps(const char * name)
{
#define FREE_MEMORY       \
	FREE_PTR(si)          \
	FREE_PTR(mat)         \
	FREE_ARRAY_PTR(clbd)  \
	FREE_ARRAY_PTR(cubd)  \
	FREE_ARRAY_PTR(ctype) \
	FREE_ARRAY_PTR(obj)   \
	FREE_ARRAY_PTR(rlbd)  \
	FREE_ARRAY_PTR(rubd)

	assert(model_);

	/** model info */
	OsiSolverInterface * si = NULL;
	CoinPackedMatrix * mat = NULL;
	double * clbd   = NULL;
	double * cubd   = NULL;
	double * obj    = NULL;
	char *   ctype  = NULL;
	double * rlbd   = NULL;
	double * rubd   = NULL;

	BGN_TRY_CATCH

	/** get DE model */
	DSP_RTN_CHECK_THROW(model_->getFullModel(mat, clbd, cubd, ctype, obj, rlbd, rubd));
	/** create solver */
	si = new OsiCbcSolverInterface();
	/** load problem */
	si->loadProblem(*mat, clbd, cubd, obj, ctype, rlbd, rubd);
	/** write mps */
	si->writeMps(name);

	/** save memory */
	FREE_MEMORY

	END_TRY_CATCH_RTN(FREE_MEMORY,;)
#undef FREE_MEMORY
}
