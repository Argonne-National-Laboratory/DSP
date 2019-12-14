/*
 * DeDriver.cpp
 *
 *  Created on: Feb 17, 2016
 *      Author: kibaekkim
 */

#include "Model/TssModel.h"
#include "Solver/Deterministic/DeDriver.h"
#include "SolverInterface/DspOsi.h"

DeDriver::DeDriver(
		DecModel *   model,  /**< model pointer */
		DspParams *  par,    /**< parameters */
		DspMessage * message /**< message pointer */):
DecSolver(model,par,message) {}

DeDriver::DeDriver(const DeDriver& rhs) :
DecSolver(rhs) {}

DeDriver::~DeDriver()
{}

/** initilize */
DSP_RTN_CODE DeDriver::init()
{
	BGN_TRY_CATCH

	primsol_.resize(model_->getFullModelNumCols());

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
				ctype[j] = 'C';
			}
		}
		if (par_->getBoolPtrParam("RELAX_INTEGRALITY")[1])
		{
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
				if (ctype[j] != 'C') {
					ctype[j] = 'C';
				}
			}
		}
	}

	switch(par_->getIntParam("SOLVER/MIP")) {
	case OsiCpx: {
#ifdef DSP_HAS_CPX
		si_ = new OsiCpxSolverInterface();
		si_->messageHandler()->setLogLevel(par_->getIntParam("LOG_LEVEL"));
		OsiCpxSolverInterface* cpx = dynamic_cast<OsiCpxSolverInterface*>(si_);
		CPXsetintparam(cpx->getEnvironmentPtr(), CPX_PARAM_THREADS, par_->getIntParam("NUM_CORES"));
#else
		throw CoinError("OsiCpx is not available.", "run", "DeDriver");
#endif
		break;
	}
	case OsiScip: {
#ifdef DSP_HAS_SCIP
		si_ = new OsiScipSolverInterface();
		si_->messageHandler()->setLogLevel(CoinMin(par_->getIntParam("LOG_LEVEL") + 2, 5));
#else
		throw CoinError("OsiScip is not available.", "run", "DeDriver");
#endif
		break;
	}
	default:
		si_ = new OsiClpSolverInterface();
		si_->messageHandler()->setLogLevel(par_->getIntParam("LOG_LEVEL"));
		break;
	}

	/** load problem */
	si_->loadProblem(*mat, clbd, cubd, obj, rlbd, rubd);
	for (int j = 0; j < mat->getNumCols(); j++)
	{
		if (ctype[j] != 'C')
			si_->setInteger(j);
	}

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
	solve();

	/** toc */
	cputime_  = CoinCpuTime() - cputime_;
	walltime_ = CoinGetTimeOfDay() - walltime_;

	/** get solutions */
	if (status_ == DSP_STAT_OPTIMAL ||
		status_ == DSP_STAT_STOPPED_TIME ||
		status_ == DSP_STAT_STOPPED_NODE ||
		status_ == DSP_STAT_STOPPED_GAP ||
	   	status_ == DSP_STAT_LIM_ITERorTIME)
	{
		/** objective bounds */
		bestprimobj_ = si_->getObjValue();
		bestdualobj_ = si_->getObjValue();

		/** solution */
		if (si_->getColSolution())
			CoinCopyN(si_->getColSolution(), si_->getNumCols(), &primsol_[0]);

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

DSP_RTN_CODE DeDriver::solve() {
	if (par_->getIntParam("SOLVER/MIP") == OsiCpx) {
#ifdef DSP_HAS_CPX
		OsiCpxSolverInterface* cpx = dynamic_cast<OsiCpxSolverInterface*>(si_);
		if (cpx) {
			cpx->branchAndBound();
		}
#endif
	} else {
		si_->initialSolve();
	}
	return DSP_RTN_OK;
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

	switch(par_->getIntParam("SOLVER/MIP")) {
	case OsiCpx:
#ifdef DSP_HAS_CPX
		si = new OsiCpxSolverInterface();
#else
		throw CoinError("OsiCpx is not available.", "writeExtMps", "DeDriver");
#endif
		break;
	case OsiScip:
#ifdef DSP_HAS_SCIP
		si = new OsiScipSolverInterface();
#else
		throw CoinError("OsiScip is not available.", "writeExtMps", "DeDriver");
#endif
		break;
	default:
		si = new OsiClpSolverInterface();
		break;
	}

	/** load problem */
	si->loadProblem(*mat, clbd, cubd, obj, ctype, rlbd, rubd);

	/** write mps */
	si->writeMps(name);

	/** save memory */
	FREE_MEMORY

	END_TRY_CATCH_RTN(FREE_MEMORY,;)
#undef FREE_MEMORY
}
