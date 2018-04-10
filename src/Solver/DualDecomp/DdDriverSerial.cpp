/*
 * DdDriverSerial.cpp
 *
 *  Created on: Jul 8, 2016
 *      Author: kibaekkim
 */

//#define DSP_DEBUG

#include "DdDriverSerial.h"
#include "Solver/DualDecomp/DdMWSerial.h"

DdDriverSerial::DdDriverSerial(DspParams* par, DecModel* model):
DdDriver(par,model)
{
}

DSP_RTN_CODE DdDriverSerial::init()
{
	BGN_TRY_CATCH

	/** create Master-Worker framework */
	mw_ = new DdMWSerial(model_, par_, message_);
	DSPdebugMessage("Created master-worker framework.\n");

	/** initialize master-worker framework */
	mw_->init();
	DSPdebugMessage("Initialized master-worker framework.\n");

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdDriverSerial::run()
{
	BGN_TRY_CATCH

	/** tic */
	cputime_  = CoinCpuTime();
	walltime_ = CoinGetTimeOfDay();

	/** run */
	DSP_RTN_CHECK_THROW(mw_->run());

	/** toc */
	cputime_  = CoinCpuTime() - cputime_;
	walltime_ = CoinGetTimeOfDay() - walltime_;

	/** retrieve master pointer */
	DdMaster * master = mw_->master_;
	if (master)
	{
		status_ = master->getStatus();
		primobj_ = master->getBestPrimalObjective();
		dualobj_ = master->getBestDualObjective();
		primsol_ = new double [model_->getFullModelNumCols()];
		dualsol_ = new double [model_->getNumCouplingRows()];
		CoinCopyN(master->getBestPrimalSolution(), model_->getFullModelNumCols(), primsol_);
		CoinCopyN(master->getBestDualSolution(), model_->getNumCouplingRows(), dualsol_);
		numNodes_ = master->getSiPtr()->getNumNodes();
		numIterations_ = master->getSiPtr()->getIterationCount();
	}
	/** nullify master pointer */
	master = NULL;

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdDriverSerial::finalize()
{
	BGN_TRY_CATCH

	/** finalize master-worker framework */
	DSP_RTN_CHECK_THROW(mw_->finalize());

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}
