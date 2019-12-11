/*
 * BdMWSerial.cpp
 *
 *  Created on: Jul 10, 2016
 *      Author: kibaekkim
 */

// #define DSP_DEBUG

#include "Solver/Benders/BdMWSerial.h"
#include "Solver/Benders/SCIPconshdlrBenders.h"
#include "SolverInterface/OsiScipSolverInterface.hpp"

BdMWSerial::BdMWSerial(
		DecModel *   model,  /**< model pointer */
		DspParams *  par,    /**< parameters */
		DspMessage * message /**< message pointer */):
BdMW(model, par, message)
{
}

BdMWSerial::~BdMWSerial() {
	// TODO Auto-generated destructor stub
}

DSP_RTN_CODE BdMWSerial::init()
{
	BGN_TRY_CATCH

	/** set parameters */
	par_->setIntPtrParamSize("ARR_PROC_IDX", model_->getNumSubproblems());
	for (int s = 0; s < model_->getNumSubproblems(); ++s)
		par_->setIntPtrParam("ARR_PROC_IDX", s, s);

	/** create master */
	master_ = new BdMaster(model_, par_, message_);

	/** initialize master */
	master_->init();

	/** create Benders worker */
	worker_ = new BdWorker(model_, par_, message_);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE BdMWSerial::run()
{
	BGN_TRY_CATCH

	/** set constraint handler */
	DSPdebugMessage("setConshdlr\n");
	master_->setConshdlr(constraintHandler());

	/** solve */
	DSPdebugMessage("solve\n");
	master_->solve();

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE BdMWSerial::finalize()
{
	BGN_TRY_CATCH

	/** finalize and free master */
	master_->finalize();
	FREE_PTR(master_);

	/** free worker */
	FREE_PTR(worker_);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

SCIPconshdlrBenders* BdMWSerial::constraintHandler()
{
	SCIPconshdlrBenders * conshdlr = NULL;

	BGN_TRY_CATCH

	/** get solver interface */
	OsiScipSolverInterface * si = dynamic_cast<OsiScipSolverInterface*>(master_->getSiPtr());

	/** Benders constraint handler */
	conshdlr = new SCIPconshdlrBenders(si->getScip(), "Benders", par_->getIntParam("BD/CUT_PRIORITY"));
	conshdlr->setDecModel(model_);
	conshdlr->setBdSub(worker_->getBdSubPtr());
	conshdlr->setOriginalVariables(si->getNumCols(),
			si->getScipVars(), par_->getIntParam("BD/NUM_CUTS_PER_ITER"));

	END_TRY_CATCH_RTN(;,NULL)

	return conshdlr;
}
