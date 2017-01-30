/*
 * DwSolverSerial.cpp
 *
 *  Created on: Dec 5, 2016
 *      Author: kibaekkim
 */

#include <DantzigWolfe/DwSolverSerial.h>
#include <DantzigWolfe/DwMaster.h>
#include <DantzigWolfe/DwMasterTr.h>
#include <DantzigWolfe/DwMasterTrLight.h>
#include <DantzigWolfe/DwWorker.h>

DwSolverSerial::DwSolverSerial(
		DecModel *   model,  /**< model pointer */
		DspParams *  par,    /**< parameters */
		DspMessage * message /**< message pointer */):
DecSolver(model, par, message),
master_(NULL),
worker_(NULL),
alps_(NULL) {}

DwSolverSerial::~DwSolverSerial() {
	FREE_PTR(master_);
	FREE_PTR(worker_);
	FREE_PTR(alps_);
}

DSP_RTN_CODE DwSolverSerial::init() {
	BGN_TRY_CATCH

	/** create worker */
	worker_ = new DwWorker(model_, par_, message_);

	/** create master */
	if (par_->getBoolParam("DW/TRUST_REGION")) {
		if (par_->getBoolParam("DW/MASTER/BRANCH_ROWS"))
			master_ = new DwMasterTr(worker_);
		else
			master_ = new DwMasterTrLight(worker_);
	} else
		master_ = new DwMaster(worker_);

	/** initialize master */
	DSP_RTN_CHECK_THROW(master_->init());

	/** create an Alps model */
	alps_ = new DspModel(master_);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

DSP_RTN_CODE DwSolverSerial::solve() {
	BGN_TRY_CATCH
	DSP_RTN_CHECK_THROW(alps_->solve());
	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

DSP_RTN_CODE DwSolverSerial::finalize() {
	BGN_TRY_CATCH

	/** finalize master */
	DSP_RTN_CHECK_THROW(master_->finalize());

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}
