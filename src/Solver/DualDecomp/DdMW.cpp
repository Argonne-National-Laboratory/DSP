/*
 * DdMW.cpp
 *
 *  Created on: Feb 9, 2016
 *      Author: kibaekkim
 */

//#define DSP_DEBUG
#include "Solver/DualDecomp/DdMW.h"

DdMW::DdMW(
		DecModel *   model,  /**< model pointer */
		DspParams *  par,    /**< parameters */
		DspMessage * message /**< message pointer */):
BaseMasterWorker(),
master_(NULL), model_(model), par_(par), message_(message),
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

	/** run master process */
	runMaster();

	/** run worker processes */
	runWorker();

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

/** initialize */
DSP_RTN_CODE DdMW::init()
{
	BGN_TRY_CATCH

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
	DSPdebug({
		DSPdebugMessage("Number of elements in xvec: %d\n", xvec->getNumElements());
		message_->printArray(xvec);});
	dup = duplicateVector(xvec, solutions);
	DSPdebugMessage("duplicateVector: %s\n", dup ? "true" : "false");

	/** free if duplicate */
	if (dup) FREE_PTR(xvec);

	return xvec;
}

void DdMW::printHeaderInfo() {
	/**
	 * PRINT DISPLAY
	 *
	 * 0123456789012345678901234567890123456789
	 * iter      curobj     primobj dualobj gap time
	 *  %4d  %+10e  %+10e  %+10e  %6.2f  %6.1f
	 */
	message_->print(1, "  %4s  %13s  %13s  %13s  %6s  %6s\n",
			"iter", "curobj", "primobj", "dualobj", "gap(%)", "times");
}

void DdMW::printIterInfo() {
	message_->print(1, " %c%4d  %+10e", itercode_, itercnt_, master_->getPrimalObjective());
	if (master_->getBestPrimalObjective() < 1.0e+20)
		message_->print(1, "  %+10e", master_->getBestPrimalObjective());
	else
		message_->print(1, "  %+13s", "Large");
	if (master_->getBestDualObjective() > -1.0e+20)
		message_->print(1, "  %+10e", master_->getBestDualObjective());
	else
		message_->print(1, "  %+13s", "Large");
	if (master_->getRelDualityGap() >= 10)
		message_->print(1, "  %6s", "Large");
	else
		message_->print(1, "  %6.2f", master_->getRelDualityGap()*100);
	message_->print(1, "  %6.1f\n", CoinGetTimeOfDay() - iterstime_);

	s_itertime_.push_back(CoinGetTimeOfDay() - iterstime_);
	s_masterobj_.push_back(master_->getPrimalObjective());
	s_bestprimobj_.push_back(master_->getBestPrimalObjective());
	s_bestdualobj_.push_back(master_->getBestDualObjective());
}

void DdMW::writeIterInfo(const char * filename)
{
	BGN_TRY_CATCH

	ofstream myfile;
	myfile.open(filename);
	myfile << "Iteration";
	myfile << ",MasterObj";
	myfile << ",BestPrimalObj";
	myfile << ",BestDualObj";
	myfile << ",Time\n";
	for (unsigned i = 0; i < s_itertime_.size(); ++i)
	{
		myfile << i;
		myfile << "," << s_masterobj_[i];
		myfile << "," << s_bestprimobj_[i];
		myfile << "," << s_bestdualobj_[i];
		myfile << "," << s_itertime_[i];
		myfile << "\n";
	}
	myfile.close();

	END_TRY_CATCH(;)
}
