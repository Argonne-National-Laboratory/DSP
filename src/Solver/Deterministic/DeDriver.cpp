/*
 * DeDriver.cpp
 *
 *  Created on: Feb 17, 2016
 *      Author: kibaekkim
 */
//#define DSP_DEBUG

#include "Model/TssModel.h"
#include "Solver/Deterministic/DeDriver.h"
#include "SolverInterface/DspOsiCpx.h"
#include "SolverInterface/DspOsiGrb.h"
#include "SolverInterface/DspOsiScip.h"

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
	FREE_PTR(qobj)		  \
	FREE_ARRAY_PTR(rlbd)  \
	FREE_ARRAY_PTR(rubd)

	assert(model_);

	/** model info */
	CoinPackedMatrix * mat = NULL;
	double * clbd   = NULL;
	double * cubd   = NULL;
	double * obj    = NULL;
	CoinPackedMatrix * qobj = NULL;
	char *   ctype  = NULL;
	double * rlbd   = NULL;
	double * rubd   = NULL;

	BGN_TRY_CATCH

	/** get DE model */
	DSP_RTN_CHECK_THROW(model_->getFullModel(mat, clbd, cubd, ctype, obj, qobj, rlbd, rubd));

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

	/** create DspOsi */
	osi_ = createDspOsi();
	if (!osi_) throw CoinError("Failed to create DspOsi", "run", "DeDriver");

	/** set display */
	osi_->setLogLevel(par_->getIntParam("DE/SOLVER/LOG_LEVEL"));

	/** set number of cores */
	osi_->setNumCores(par_->getIntParam("NUM_CORES"));

	/** load problem */
	osi_->si_->loadProblem(*mat, clbd, cubd, obj, rlbd, rubd);
	//PRINT_ARRAY_MSG(qobj->getNumElements(), qobj->getElements(), "elements in qobj");
	if (qobj != NULL){
		osi_->loadQuadraticObjective(*qobj);
	}
	
	for (int j = 0; j < mat->getNumCols(); j++)
	{
		if (ctype[j] != 'C')
			osi_->si_->setInteger(j);
	}

	/** set optimality gap tolerance */
	osi_->setRelMipGap(par_->getDblParam("DE/GAPTOL"));

	/** time limit */
	double time_limit = par_->getDblParam("DE/WALL_LIM");
	osi_->setTimeLimit(CoinMin(time_remains_,time_limit));

	/** set node limit */
	osi_->setNodeLimit(par_->getIntParam("NODE_LIM"));

	/** tic */
	cputime_  = CoinCpuTime();
	walltime_ = CoinGetTimeOfDay();

	/** solve */
	solve();

	/** toc */
	cputime_  = CoinCpuTime() - cputime_;
	walltime_ = CoinGetTimeOfDay() - walltime_;

	/** get solutions */
	status_ = osi_->status();
	if (status_ == DSP_STAT_OPTIMAL ||
		status_ == DSP_STAT_STOPPED_TIME ||
		status_ == DSP_STAT_STOPPED_NODE ||
		status_ == DSP_STAT_STOPPED_GAP ||
	   	status_ == DSP_STAT_LIM_ITERorTIME)
	{
		/** objective bounds */
		bestprimobj_ = osi_->getPrimObjValue();
		bestdualobj_ = osi_->getDualObjValue();

		/** solution */
		if (osi_->si_->getColSolution()) {
			CoinCopyN(osi_->si_->getColSolution(), osi_->si_->getNumCols(), &primsol_[0]);
			bestprimsol_ = primsol_;
		}

		/** statistics */
		numIterations_ = osi_->si_->getIterationCount();
		numNodes_ = osi_->getNumNodes();
	}
	// osi_->si_->writeMps("dsp");

	/** save memory */
	FREE_MEMORY

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	return DSP_RTN_OK;

#undef FREE_MEMORY
}

DSP_RTN_CODE DeDriver::solve() {
	osi_->solve();
	return DSP_RTN_OK;
}

DSP_RTN_CODE DeDriver::finalize() {
	return DSP_RTN_OK;
}

void DeDriver::writeExtMps(const char * name)
{
#define FREE_MEMORY       \
	FREE_PTR(osi)         \
	FREE_PTR(mat)         \
	FREE_ARRAY_PTR(clbd)  \
	FREE_ARRAY_PTR(cubd)  \
	FREE_ARRAY_PTR(ctype) \
	FREE_ARRAY_PTR(obj)   \
	FREE_ARRAY_PTR(rlbd)  \
	FREE_ARRAY_PTR(rubd)

	assert(model_);

	/** model info */
	DspOsi * osi = NULL;
	CoinPackedMatrix * mat = NULL;
	double * clbd   = NULL;
	double * cubd   = NULL;
	double * obj    = NULL;
	CoinPackedMatrix * qobj = NULL;
	char *   ctype  = NULL;
	double * rlbd   = NULL;
	double * rubd   = NULL;
	
	BGN_TRY_CATCH

	/** get DE model */
	DSP_RTN_CHECK_THROW(model_->getFullModel(mat, clbd, cubd, ctype, obj, qobj, rlbd, rubd));

	/** create DspOsi */
	osi = createDspOsi();
	if (!osi) throw CoinError("Failed to create DspOsi", "writeExtMps", "DeDriver");

	/** load problem */
	osi->si_->loadProblem(*mat, clbd, cubd, obj, rlbd, rubd);
	if (qobj!=NULL){
		osi->loadQuadraticObjective(*qobj);
	}
	
	for (int j = 0; j < mat->getNumCols(); j++)
	{
		if (ctype[j] != 'C')
			osi->si_->setInteger(j);
	}

	/** write mps */
	osi->writeMps(name);

	/** save memory */
	FREE_MEMORY

	END_TRY_CATCH_RTN(FREE_MEMORY,;)
#undef FREE_MEMORY
}

DspOsi * DeDriver::createDspOsi() {
	DspOsi * osi = NULL;
	BGN_TRY_CATCH

	switch(par_->getIntParam("DE/SOLVER")) {
	case OsiCpx:
#ifdef DSP_HAS_CPX
		osi = new DspOsiCpx();
#else
		throw CoinError("Cplex is not available.", "createDspOsi", "DeDriver");
#endif
		break;
	case OsiGrb:
#ifdef DSP_HAS_GRB
		osi = new DspOsiGrb();
#else
		throw CoinError("Gurobi is not available.", "createDspOsi", "DeDriver");
#endif
		break;
	case OsiScip:
#ifdef DSP_HAS_SCIP
		osi = new DspOsiScip();
#else
		throw CoinError("Scip is not available.", "createDspOsi", "DeDriver");
#endif
		break;
	default:
		throw CoinError("Invalid paramter value", "createDspOsi", "DeDriver");
		break;
	}

	END_TRY_CATCH_RTN(;,osi)
	return osi;
}
