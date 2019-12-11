/*
 * BdMaster.cpp
 *
 *  Created on: Feb 18, 2016
 *      Author: kibaekkim
 */

//#define DSP_DEBUG

#include "Model/TssModel.h"
#include "Solver/Benders/BdMaster.h"
#include "SolverInterface/OsiScipSolverInterface.hpp"
#include "Solver/Benders/SCIPconshdlrBendersWorker.h"

BdMaster::BdMaster(
			DecModel *   model,   /**< model pointer */
			DspParams *  par,     /**< parameter pointer */
			DspMessage * message /**< message pointer */) :
DecSolver(model, par, message),
worker_(NULL),
naux_(1),
obj_aux_(NULL),
clbd_aux_(NULL),
cubd_aux_(NULL) {
	obj_aux_ = new double[naux_];
	clbd_aux_ = new double[naux_];
	cubd_aux_ = new double[naux_];
	obj_aux_[0] = 1.0;
	clbd_aux_[0] = -COIN_DBL_MAX;
	cubd_aux_[0] = +COIN_DBL_MAX;
	tic_ = CoinGetTimeOfDay();
}

BdMaster::BdMaster(const BdMaster& rhs) :
DecSolver(rhs),
naux_(rhs.naux_) {
	obj_aux_ = new double[naux_];
	clbd_aux_ = new double[naux_];
	cubd_aux_ = new double[naux_];
	for (int j = 0; j < naux_; ++j) {
		obj_aux_[j] = rhs.obj_aux_[j];
		clbd_aux_[j] = rhs.clbd_aux_[j];
		cubd_aux_[j] = rhs.cubd_aux_[j];
	}
}

BdMaster::~BdMaster() {
	FREE_ARRAY_PTR(obj_aux_);
	FREE_ARRAY_PTR(clbd_aux_);
	FREE_ARRAY_PTR(cubd_aux_);
}

DSP_RTN_CODE BdMaster::init() {
	BGN_TRY_CATCH

	/** create problem */
	createProblem();

	/** set node limit */
	si_->setNodeLimit(par_->getIntParam("NODE_LIM"));

	/** set print level */
	si_->messageHandler()->setLogLevel(CoinMin(par_->getIntParam("LOG_LEVEL") + 2, 5));

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE BdMaster::solve() {

	BGN_TRY_CATCH

	DSPdebugMessage("Start solving...\n");

	/** set time limit */
	si_->setTimeLimit(CoinMax(0.01, CoinMin(1.0e+20,time_remains_)));
	tic_ = CoinGetTimeOfDay();

	/** solve */
	si_->initialSolve();

	/** solver status */
	status_ = getStatus();
	DSPdebugMessage("Benders status %d\n", status_);
	switch(status_)
	{
	case DSP_STAT_OPTIMAL:
	case DSP_STAT_LIM_ITERorTIME:
	case DSP_STAT_STOPPED_GAP:
	case DSP_STAT_STOPPED_NODE:
	case DSP_STAT_STOPPED_TIME:
	{
		/** get solution */
		if (si_->getColSolution() != NULL)
			CoinCopyN(si_->getColSolution(), si_->getNumCols(), &primsol_[0]);
		/** primal objective value */
		primobj_ = si_->getObjValue();
		/** dual objective value */
		dualobj_ = si_->getBestDualBound();
		break;
	}
	default:
		message_->print(0, "Warning: master solution status is %d\n", status_);
		break;
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE BdMaster::setSolutions(Solutions initsols)
{
	BGN_TRY_CATCH

	if (si_->getNumIntegers() > 0)
	{
		OsiScipSolverInterface * SiScip = dynamic_cast<OsiScipSolverInterface*>(si_);
		for (unsigned i = 0; i < initsols.size(); ++i)
		{
			double * solution = initsols[i]->denseVector(SiScip->getNumCols());
			SiScip->setColSolution(solution);
			FREE_ARRAY_PTR(solution);
		}
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE BdMaster::setConshdlr(SCIPconshdlrBenders* conshdlr)
{
	BGN_TRY_CATCH

	/** retrieve solver interface for SCIP */
	OsiScipSolverInterface * scip = dynamic_cast<OsiScipSolverInterface*>(si_);

	/** include constraint handler */
	SCIP_CALL_ABORT(SCIPincludeObjConshdlr(scip->getScip(), conshdlr, false));

	/** create constraint */
	SCIP_CONS * cons = NULL;
	SCIP_CALL_ABORT(SCIPcreateConsBenders(scip->getScip(), &cons, "Benders"));
	SCIP_CALL_ABORT(SCIPaddCons(scip->getScip(), cons));
	SCIP_CALL_ABORT(SCIPreleaseCons(scip->getScip(), &cons));
	DSPdebugMessage("Added constraint handler %p\n", conshdlr);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE BdMaster::createProblem() {
#define FREE_MEMORY       \
	FREE_PTR(mat)         \
	FREE_ARRAY_PTR(clbd)  \
	FREE_ARRAY_PTR(cubd)  \
	FREE_ARRAY_PTR(ctype) \
	FREE_ARRAY_PTR(obj)   \
	FREE_ARRAY_PTR(rlbd)  \
	FREE_ARRAY_PTR(rubd)

	assert(model_);

	if (naux_ <= 0 || !obj_aux_ || !clbd_aux_ || !cubd_aux_) {
		printf("Warning: Auxiliary column information is required.\n");
		return DSP_RTN_ERR;
	}

	/** master problem */
	CoinPackedMatrix * mat = NULL;
	double * clbd   = NULL;
	double * cubd   = NULL;
	double * obj    = NULL;
	char *   ctype  = NULL;
	double * rlbd   = NULL;
	double * rubd   = NULL;

	BGN_TRY_CATCH

	/** number of columns */
	int ncols = model_->getNumCouplingCols() + naux_;
	
	/** number of integer variables in the core */
	int nIntegers = model_->getNumIntegers();
	DSPdebugMessage("ncols %d, nIntegers %d\n", ncols, nIntegers);

	/** decompose model */
	DSP_RTN_CHECK_THROW(model_->decompose(
			0, NULL, naux_, clbd_aux_, cubd_aux_, obj_aux_,
			mat, clbd, cubd, ctype, obj, rlbd, rubd));
	DSPdebugMessage("Decomposed the model.\n");

	/** convert column types */
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
			DSPdebugMessage("Relax first-stage integrality.\n");
			for (int j = 0; j < tssModel->getNumCols(0); ++j)
			{
				if (ctype[j] != 'C')
					nIntegers--;
				ctype[j] = 'C';
			}
		}
	}
	else
	{
		if (par_->getBoolPtrParam("RELAX_INTEGRALITY")[0])
		{
			for (int j = 0; j < mat->getNumCols(); j++)
			{
				if (ctype[j] != 'C')
					nIntegers--;
				ctype[j] = 'C';
			}
		}
	}

//	if (nIntegers > 0)
//	{
		assert(si_==NULL);
		si_ = new OsiScipSolverInterface();
		si_->messageHandler()->setLogLevel(CoinMin(par_->getIntParam("LOG_LEVEL"), 5));
		DSPdebugMessage("Successfully created SCIP interface \n");
//	}
//	else
//	{
//		si_ = new SolverInterfaceClp(par_);
//		si_->setPrintLevel(CoinMax(par_->getIntParam("LOG_LEVEL") - 1, 0));
//		DSPdebugMessage("Successfully created Soplex interface \n");
//	}

	/** load problem data */
	si_->loadProblem(*mat, clbd, cubd, obj, rlbd, rubd);
	for (int j = 0; j < mat->getNumCols(); ++j) {
		if (ctype[j] != 'C')
			si_->setInteger(j);
	}
	DSPdebugMessage("Successfully load the problem.\n");

	/** allocate memory for primal solution */
	primsol_.resize(si_->getNumCols());

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	/** save memory */
	FREE_MEMORY

	return DSP_RTN_OK;

#undef FREE_MEMORY
}

DSP_RTN_CODE BdMaster::setObjectiveBounds(double upper, double lower) {
#define FREE_MEMORY         \
	FREE_ARRAY_PTR(auxind)  \
	FREE_ARRAY_PTR(auxcoef)

	/** for recourse lower bound */
	int * auxind     = NULL;
	double * auxcoef = NULL;

	BGN_TRY_CATCH

	/** number of columns */
	int ncols = si_->getNumCols();

	/** allocate memory */
	auxind = new int [ncols];
	auxcoef = new double [ncols];

	/** update bounds */
	primobj_ = upper;
	dualobj_ = lower;

	for (int j = 0; j < ncols; ++j) {
		auxind[j] = j;
		auxcoef[j] = si_->getObjCoefficients()[j];
	}
	si_->addRow(ncols, auxind, auxcoef, dualobj_, primobj_);

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

DSP_RTN_CODE BdMaster::setAuxVarData(
		int size,
		double* obj,
		double* clbd,
		double* cubd) {
	BGN_TRY_CATCH

	/** free memory (just in case) */
	FREE_ARRAY_PTR(obj_aux_)
	FREE_ARRAY_PTR(clbd_aux_)
	FREE_ARRAY_PTR(cubd_aux_)
	/** allocate memory */
	naux_ = size;
	obj_aux_ = new double[naux_];
	clbd_aux_ = new double[naux_];
	cubd_aux_ = new double[naux_];
	/** copy data */
	CoinCopyN(obj, naux_, obj_aux_);
	CoinCopyN(clbd, naux_, clbd_aux_);
	CoinCopyN(cubd, naux_, cubd_aux_);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}
