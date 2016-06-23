/*
 * BdWorker.cpp
 *
 *  Created on: Feb 18, 2016
 *      Author: kibaekkim
 */

#define DSP_DEBUG

#include "Solver/Benders/BdWorker.h"

BdWorker::BdWorker(DspParams * par, DecModel * model, DspMessage * message):
	DecSolver(par, model, message),
	bdsub_(NULL),
	parProcIdxSize_(0), parProcIdx_(NULL) {}

BdWorker::~BdWorker()
{
	FREE_PTR(bdsub_);
}

DSP_RTN_CODE BdWorker::init()
{
	BGN_TRY_CATCH

	/** parameters */
	parProcIdxSize_ = par_->getIntPtrParamSize("ARR_PROC_IDX");
	parProcIdx_     = par_->getIntPtrParam("ARR_PROC_IDX");
	DSPdebugMessage("parProcIdxSize_ %d\n", parProcIdxSize_);

	/** create problem */
	createProblem();

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE BdWorker::solve()
{
	BGN_TRY_CATCH

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE BdWorker::createProblem()
{
	/** subproblems */
	TssModel * tssmodel = NULL;

	BGN_TRY_CATCH

	tssmodel = dynamic_cast<TssModel*>(model_);
	if (!tssmodel) throw "Invalid model type cast";

	bdsub_ = new BdSub(par_);
	DSP_RTN_CHECK_THROW(bdsub_->setSubIndices(parProcIdxSize_, parProcIdx_), "setSubIndices", "BdSub");
	DSP_RTN_CHECK_THROW(bdsub_->loadProblem(tssmodel), "loadProblem", "BdSub");

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}
