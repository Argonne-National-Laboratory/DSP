/*
 * DeDriver.cpp
 *
 *  Created on: Feb 17, 2016
 *      Author: kibaekkim
 */

#include "Model/TssModel.h"
#include "Solver/Deterministic/DeDriver.h"
#include "SolverInterface/SolverInterfaceClp.h"

#ifndef NO_CPX
#include "SolverInterface/SolverInterfaceCpx.h"
#endif

#ifndef NO_SCIP
#include "SolverInterface/SolverInterfaceScip.h"
#endif

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

	if (nIntegers > 0)
	{
		switch(par_->getIntParam("SOLVER/MIP")) {
		case CPLEX: {
#ifndef NO_CPX
			si_ = new SolverInterfaceCpx(par_);
			si_->setPrintLevel(par_->getIntParam("LOG_LEVEL"));
			SolverInterfaceCpx* osicpx = dynamic_cast<SolverInterfaceCpx*>(si_);
			OsiCpxSolverInterface* cpx = dynamic_cast<OsiCpxSolverInterface*>(osicpx->getOSI());
			CPXsetintparam(cpx->getEnvironmentPtr(), CPX_PARAM_THREADS, par_->getIntParam("NUM_CORES"));
			break;
#endif
		}
		case SCIP: {
#ifndef NO_SCIP
			si_ = new SolverInterfaceScip(par_);
			si_->setPrintLevel(CoinMin(par_->getIntParam("LOG_LEVEL") + 2, 5));
			break;
#endif
		}
		default:
			break;
		}
	}
	else
	{
		si_ = new SolverInterfaceClp(par_);
		/** print level */
		si_->setPrintLevel(par_->getIntParam("LOG_LEVEL"));
	}

	/** load problem */
	si_->loadProblem(mat, clbd, cubd, obj, ctype, rlbd, rubd, "DeDriver");

	/** time limit */
	double time_limit = CoinMin(
			par_->getDblParam("DE/WALL_LIM"),
			par_->getDblParam("MIP/TIME_LIM"));
	si_->setTimeLimit(time_limit);

	/** set node limit */
	si_->setNodeLimit(par_->getIntParam("NODE_LIM"));

	/** tic */
	cputime_  = CoinCpuTime();
	walltime_ = CoinGetTimeOfDay();

	/** solve */
	si_->solve();

	/** toc */
	cputime_  = CoinCpuTime() - cputime_;
	walltime_ = CoinGetTimeOfDay() - walltime_;

	/** solution status */
	status_ = si_->getStatus();

	/** get solutions */
	if (status_ == DSP_STAT_OPTIMAL ||
		status_ == DSP_STAT_STOPPED_TIME ||
		status_ == DSP_STAT_STOPPED_NODE ||
		status_ == DSP_STAT_STOPPED_GAP)
	{
		/** objective bounds */
		primobj_ = si_->getPrimalBound();
		dualobj_ = si_->getDualBound();

		/** solution */
		if (si_->getSolution())
			CoinCopyN(si_->getSolution(), si_->getNumCols(), primsol_);

		/** statistics */
		numIterations_ = si_->getIterationCount();
		numNodes_ = si_->getNumNodes();
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
	SolverInterface * si = NULL;
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

	switch(par_->getIntParam("SOLVER/MIP")) {
	case CPLEX:
#ifndef NO_CPX
		si = new SolverInterfaceCpx(par_);
		break;
#endif
	case SCIP:
#ifndef NO_SCIP
		si = new SolverInterfaceScip(par_);
		break;
#endif
	default:
		si = new SolverInterfaceClp(par_);
		break;
	}

	/** load problem */
	si->loadProblem(mat, clbd, cubd, obj, ctype, rlbd, rubd, "writeMps");

	/** write mps */
	si->writeMps(name);

	/** save memory */
	FREE_MEMORY

	END_TRY_CATCH_RTN(FREE_MEMORY,;)
#undef FREE_MEMORY
}
