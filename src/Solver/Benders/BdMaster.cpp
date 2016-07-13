/*
 * BdMaster.cpp
 *
 *  Created on: Feb 18, 2016
 *      Author: kibaekkim
 */

//#define DSP_DEBUG
#include "Solver/Benders/BdMaster.h"
#include "SolverInterface/SolverInterfaceScip.h"
#include "SolverInterface/SolverInterfaceClp.h"
#include "Solver/Benders/SCIPconshdlrBendersWorker.h"

BdMaster::BdMaster(
		DspParams* par,
		DecModel* model,
		DspMessage * message):
DecSolver(par, model, message),
si_(NULL),
naux_(1),
obj_aux_(NULL),
clbd_aux_(NULL),
cubd_aux_(NULL)
{
	obj_aux_  = new double [naux_];
	clbd_aux_ = new double [naux_];
	cubd_aux_ = new double [naux_];
	obj_aux_[0] = 1.0;
	clbd_aux_[0] = -COIN_DBL_MAX;
	cubd_aux_[0] = +COIN_DBL_MAX;
	tic_ = CoinGetTimeOfDay();
}

BdMaster::~BdMaster()
{
	FREE_PTR(si_);
	FREE_ARRAY_PTR(obj_aux_);
	FREE_ARRAY_PTR(clbd_aux_);
	FREE_ARRAY_PTR(cubd_aux_);
}

DSP_RTN_CODE BdMaster::init()
{
	BGN_TRY_CATCH

	/** create problem */
	createProblem();

	/** set node limit */
	si_->setNodeLimit(par_->getIntParam("NODE_LIM"));

	/** set print level */
	si_->setPrintLevel(CoinMin(par_->getIntParam("LOG_LEVEL") + 2, 5));

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE BdMaster::solve()
{
	BGN_TRY_CATCH

	DSPdebugMessage("Start solving...\n");

	/** set time limit */
	time_remains_ = par_->getDblParam("SCIP/TIME_LIM");
	time_remains_ -= CoinGetTimeOfDay() - tic_;
	si_->setTimeLimit(time_remains_);
	tic_ = CoinGetTimeOfDay();

	/** solve */
	si_->solve();

	/** solver status */
	status_ = si_->getStatus();
	switch(status_)
	{
	case DSP_STAT_OPTIMAL:
	case DSP_STAT_LIM_ITERorTIME:
	case DSP_STAT_STOPPED_GAP:
	case DSP_STAT_STOPPED_NODE:
	case DSP_STAT_STOPPED_TIME:
		/** get solution */
		CoinCopyN(si_->getSolution(), si_->getNumCols(), primsol_);
		/** primal objective value */
		primobj_ = si_->getPrimalBound();
		/** dual objective value */
		dualobj_ = si_->getDualBound();
		break;
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
		SolverInterfaceScip * SiScip = dynamic_cast<SolverInterfaceScip*>(si_);
		for (unsigned i = 0; i < initsols.size(); ++i)
		{
			double * solution = initsols[i]->denseVector(SiScip->getNumCols());
			SiScip->setSolution(solution);
			FREE_ARRAY_PTR(solution);
		}
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE BdMaster::setBranchingPriority(
		int   size,      /**< size of array */
		int * priorities /**< branch priority */)
{
	BGN_TRY_CATCH

	if (si_->getNumIntegers() > 0)
	{
		SolverInterfaceScip * SiScip = dynamic_cast<SolverInterfaceScip*>(si_);
		SiScip->setBranchPriorities(size, priorities);
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE BdMaster::setConshdlr(SCIPconshdlrBenders* conshdlr)
{
	BGN_TRY_CATCH

	/** retrieve solver interface for SCIP */
	SolverInterfaceScip * SiScip = dynamic_cast<SolverInterfaceScip*>(si_);

	/** add constraint handler */
	SiScip->addConstraintHandler(conshdlr, false);
	DSPdebugMessage("Added constraint handler %p\n", conshdlr);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE BdMaster::createProblem()
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

	if (naux_ <= 0 || !obj_aux_ || !clbd_aux_ || !cubd_aux_)
	{
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

	int ncols = model_->getNumSubproblemCouplingCols(0) + naux_;

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
		si_ = new SolverInterfaceScip(par_);
		si_->setPrintLevel(CoinMin(par_->getIntParam("LOG_LEVEL"), 5));
		DSPdebugMessage("Successfully created SCIP interface \n");
//	}
//	else
//	{
//		si_ = new SolverInterfaceClp(par_);
//		si_->setPrintLevel(CoinMax(par_->getIntParam("LOG_LEVEL") - 1, 0));
//		DSPdebugMessage("Successfully created Soplex interface \n");
//	}

	/** load problem data */
	si_->loadProblem(mat, clbd, cubd, obj, ctype, rlbd, rubd, "BdMaster");
	DSPdebugMessage("Successfully load the problem.\n");

	/** allocate memory for primal solution */
	primsol_ = new double [si_->getNumCols()];

	/** save memory */
	FREE_MEMORY

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	return DSP_RTN_OK;

#undef FREE_MEMORY
}

DSP_RTN_CODE BdMaster::setObjectiveBounds(double upper, double lower)
{
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

	for (int j = 0; j < ncols; ++j)
	{
		auxind[j] = j;
		auxcoef[j] = si_->getObjCoef()[j];
	}
	si_->addRow(ncols, auxind, auxcoef, dualobj_, primobj_);

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

DSP_RTN_CODE BdMaster::setAuxVarData(int size, double* obj, double* clbd, double* cubd)
{
	BGN_TRY_CATCH

	FREE_ARRAY_PTR(obj_aux_)
	FREE_ARRAY_PTR(clbd_aux_)
	FREE_ARRAY_PTR(cubd_aux_)

	naux_ = size;
	obj_aux_  = new double [naux_];
	clbd_aux_ = new double [naux_];
	cubd_aux_ = new double [naux_];

	CoinCopyN(obj, naux_, obj_aux_);
	CoinCopyN(clbd, naux_, clbd_aux_);
	CoinCopyN(cubd, naux_, cubd_aux_);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}
