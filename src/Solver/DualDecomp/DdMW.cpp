/*
 * DdMW.cpp
 *
 *  Created on: Feb 9, 2016
 *      Author: kibaekkim
 */

//#define DSP_DEBUG
#include "Solver/DualDecomp/DdMW.h"

DdMW::DdMW(
		DdMaster *        master, /**< master problem */
		vector<DdWorker*> worker  /**< worker for finding lower bounds */):
BaseMasterWorker(), master_(master), worker_(worker),
model_(NULL), par_(NULL), message_(NULL),
parFeasCuts_(-1), parOptCuts_(-1), parEvalUb_(-1),
itercode_(' '), itercnt_(0), iterstime_(0.0)
{
	cutsToAdd_ = new OsiCuts;
}

DdMW::~DdMW() {
	FREE_PTR(cutsToAdd_);
	model_   = NULL;
	par_     = NULL;
	message_ = NULL;
}

DSP_RTN_CODE DdMW::run()
{
	BGN_TRY_CATCH

	/** initialize */
	init();

	/** run master process */
	runMaster();

	/** run worker processes */
	runWorker();

	/** finalize */
	finalize();

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

/** initialize */
DSP_RTN_CODE DdMW::init() {
	BGN_TRY_CATCH

	/** retrieve model, parameter, and message objects */
	if (master_)
	{
		model_   = master_->getModelPtr();
		par_     = master_->getParPtr();
		message_ = master_->getMessagePtr();
	}
	if (worker_.size() > 0)
	{
		model_   = worker_[0]->getModelPtr();
		par_     = worker_[0]->getParPtr();
		message_ = worker_[0]->getMessagePtr();
	}

	parFeasCuts_ = par_->getIntParam("DD/FEAS_CUTS");
	parOptCuts_  = par_->getIntParam("DD/OPT_CUTS");
	parEvalUb_   = par_->getIntParam("DD/EVAL_UB");

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMW::finalize() {
	BGN_TRY_CATCH

	/** delete local solutions */
	for (unsigned i = 0; i < ubSolutions_.size(); ++i)
		FREE_PTR(ubSolutions_[i]);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

CoinPackedVector * DdMW::duplicateSolution(
		int size,           /**< size of array */
		const double * x,   /**< current solution */
		Solutions solutions /**< solution pool to check duplication */)
{
	assert(x);
	bool dup = false;

	CoinPackedVector * xvec = new CoinPackedVector;
	for (int i = 0; i < size; ++i)
	{
		if (fabs(x[i]) > 1.0e-8)
			xvec->insert(i, x[i]);
	}
#ifdef DSP_DEBUG
	DSPdebugMessage("Number of elements in xvec: %d\n", xvec->getNumElements());
	message_->printArray(xvec);
#endif
	dup = duplicateVector(xvec, solutions);
	DSPdebugMessage("duplicateVector: %s\n", dup ? "true" : "false");

	/** free if duplicate */
	if (dup) FREE_PTR(xvec);

	return xvec;
}

void DdMW::printIterInfo() {
	message_->print(1, " %c%4d  %+10e", itercode_, itercnt_, master_->getPrimalObjective());
	if (master_->getBestPrimalObjective() < 1.0e+20)
		message_->print(1, "  %+10e", master_->getBestPrimalObjective());
	else
		message_->print(1, "  %+10s", "inf");
	message_->print(1, "  %+10e  %6.2f  %6.1f\n",
			master_->getBestDualObjective(), master_->getRelDualityGap()*100, CoinGetTimeOfDay() - iterstime_);
}
