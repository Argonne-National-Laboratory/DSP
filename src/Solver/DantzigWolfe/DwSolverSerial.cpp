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
DecSolver(model, par, message), master_(NULL), worker_(NULL) {
}

DwSolverSerial::~DwSolverSerial() {
	FREE_PTR(master_);
	FREE_PTR(worker_);
}

DSP_RTN_CODE DwSolverSerial::init() {
	BGN_TRY_CATCH
	worker_ = new DwWorker(model_, par_, message_);
	if (par_->getBoolParam("DW/TRUST_REGION")) {
		if (par_->getBoolParam("DW/MASTER/BRANCH_ROWS"))
			master_ = new DwMasterTr(worker_);
		else
			master_ = new DwMasterTrLight(worker_);
	} else
		master_ = new DwMaster(worker_);
	DSP_RTN_CHECK_THROW(master_->init());
	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

DSP_RTN_CODE DwSolverSerial::solve() {
	BGN_TRY_CATCH
	DSP_RTN_CHECK_THROW(master_->solve());
	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

DSP_RTN_CODE DwSolverSerial::finalize() {
	BGN_TRY_CATCH
	DSP_RTN_CHECK_THROW(master_->finalize());
	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

bool DwSolverSerial::chooseBranchingObjects(
		DspBranch*& branchingUp, /**< [out] branching-up object */
		DspBranch*& branchingDn  /**< [out] branching-down object */) {
	return master_->chooseBranchingObjects(branchingUp, branchingDn);
}

void DwSolverSerial::setBranchingObjects(const DspBranch* branchobj) {
	master_->setBranchingObjects(branchobj);
}
