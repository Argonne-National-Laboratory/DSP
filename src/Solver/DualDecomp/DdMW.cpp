/*
 * DdMW.cpp
 *
 *  Created on: Feb 9, 2016
 *      Author: kibaekkim
 */

//#define DSP_DEBUG
#include "Utility/DspUtility.h"
#include "Solver/DualDecomp/DdMW.h"

DdMW::DdMW(
		DecModel *   model,  /**< model pointer */
		DspParams *  par,    /**< parameters */
		DspMessage * message /**< message pointer */):
BaseMasterWorker(),
master_(NULL), model_(model), par_(par), message_(message),
parFeasCuts_(-1), parOptCuts_(-1), parEvalUb_(-1), parTimeLimit_(0),
itercode_(' '), itercnt_(0), iterstime_(0.0)
{
	cutsToAdd_ = new OsiCuts;
}

DdMW::~DdMW() {
	FREE_PTR(cutsToAdd_);
	model_   = NULL;
	par_     = NULL;
	message_ = NULL;

	/** free master */
	FREE_PTR(master_);

	/** free workers */
	for (unsigned i = 0; i < worker_.size(); ++i)
		FREE_PTR(worker_[i]);
	worker_.clear();
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

	parFeasCuts_  = par_->getIntParam("DD/FEAS_CUTS");
	parOptCuts_   = par_->getIntParam("DD/OPT_CUTS");
	parEvalUb_    = par_->getIntParam("DD/EVAL_UB");
	parTimeLimit_ = par_->getDblParam("DD/WALL_LIM");

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

DSP_RTN_CODE DdMW::storeCouplingSolutions(Solutions& stored)
{
	BGN_TRY_CATCH

		/** maximum number of solutions to evaluate */
		int max_stores = par_->getIntParam("DD/MAX_EVAL_UB");

	/** store solutions to distribute */
	for (int s = 0; s < model_->getNumSubproblems(); ++s)
	{
		int nx = model_->getNumSubproblemCouplingCols(s);

		DSPdebugMessage2("ubSolutions_ %lu\n", ubSolutions_.size());
		DSPdebugMessage2("solution[%d] nx %d:\n", s, nx);
		DSPdebug2(message_->printArray(nx, master_->subsolution_[s]));

		CoinPackedVector * x = duplicateSolution(
				nx, master_->subsolution_[s], ubSolutions_);
		if (x != NULL)
		{
			DSPdebugMessage2("Coupling solution:\n");
			DSPdebug2(DspMessage::printArray(nx, master_->subsolution_[s]));
			/** store solution */
			ubSolutions_.push_back(x);
			stored.push_back(x);
			/** count */
			max_stores--;
		}
		if (max_stores == 0) break;
	}

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
	message_->print(1, "\nDUAL DECOMPOSITION ITERATION INFORMATION:\n");
	message_->print(1, "* master   = objective function value of the master problem.\n");
	message_->print(1, "* primobj  = best primal objective function value.\n");
	message_->print(1, "* dualobj  = best dual objective function value.\n");
	message_->print(1, "* a.gap(%) = Approximate gap between master and dualobj.\n");
	message_->print(1, "* d.gap(%) = Duality gap between primobj and dualobj.\n");
	message_->print(1, "* times    = wall clock time in seconds.\n\n");
	message_->print(1, "  %4s  %13s  %13s  %13s  %8s  %8s  %6s\n",
			"iter", "master", "primobj", "dualobj", "a.gap(%)", "d.gap(%)", "time");
}

void DdMW::printIterInfo()
{
	double primobj = master_->getPrimalObjective();
	double bestprimobj = master_->getBestPrimalObjective();
	double bestdualobj = master_->getBestDualObjective();
	double approxgap = master_->getRelApproxGap();
	double dualitygap = master_->getRelDualityGap();

	message_->print(1, " %c%4d", itercode_, itercnt_);
	if (primobj < 1.0e+20)
		message_->print(1, "  %+10e", primobj);
	else
		message_->print(1, "  %+13s", "Large");
	if (bestprimobj < 1.0e+20)
		message_->print(1, "  %+10e", bestprimobj);
	else
		message_->print(1, "  %+13s", "Large");
	if (bestdualobj > -1.0e+20)
		message_->print(1, "  %+10e", bestdualobj);
	else
		message_->print(1, "  %+13s", "Large");
	message_->print(1, "  %8.2f", approxgap*100);
	if (dualitygap < 10.0)
		message_->print(1, "  %8.2f", dualitygap*100);
	else
		message_->print(1, "  %8s", "Large");
	message_->print(1, "  %6.1f\n", CoinGetTimeOfDay() - iterstime_);

	s_itertime_.push_back(CoinGetTimeOfDay() - iterstime_);
	s_masterobj_.push_back(primobj);
	s_bestprimobj_.push_back(bestprimobj);
	s_bestdualobj_.push_back(bestdualobj);
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
