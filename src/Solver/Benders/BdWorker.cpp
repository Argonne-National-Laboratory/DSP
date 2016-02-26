/*
 * BdWorker.cpp
 *
 *  Created on: Feb 18, 2016
 *      Author: kibaekkim
 */

#include "Solver/Benders/BdWorker.h"

BdWorker::BdWorker(DspParams * par, DecModel * model, StoMessage * message):
	DecSolver(par, model, message),
	bdsub_(NULL),
	parProcIdxSize_(0), parProcIdx_(NULL) {}

BdWorker::~BdWorker()
{
	FREE_PTR(bdsub_);
}

STO_RTN_CODE BdWorker::init()
{
	BGN_TRY_CATCH

	/** parameters */
	parProcIdxSize_ = par_->getIntPtrParamSize("ARR_PROC_IDX");
	parProcIdx_     = par_->getIntPtrParam("ARR_PROC_IDX");

	/** create problem */
	createProblem();

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
}

STO_RTN_CODE BdWorker::solve()
{
	BGN_TRY_CATCH

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
}

STO_RTN_CODE BdWorker::createProblem()
{
	/** subproblems */
	TssModel * tssmodel = NULL;

	BGN_TRY_CATCH

	tssmodel = dynamic_cast<TssModel*>(model_);
	if (!tssmodel) throw "Invalid model type cast";

	bdsub_ = new BdSub(par_);
	STO_RTN_CHECK_THROW(bdsub_->setSubIndices(parProcIdxSize_, parProcIdx_), "setSubIndices", "BdSub");
	STO_RTN_CHECK_THROW(bdsub_->loadProblem(tssmodel), "loadProblem", "BdSub");

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
}
