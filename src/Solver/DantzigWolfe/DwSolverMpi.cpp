/*
 * DwSolverMpi.cpp
 *
 *  Created on: Dec 5, 2016
 *      Author: kibaekkim
 */

#define DSP_DEBUG
#include <DantzigWolfe/DwSolverMpi.h>
#include <DantzigWolfe/DwMasterTr.h>
#include <DantzigWolfe/DwWorkerMpi.h>

DwSolverMpi::DwSolverMpi(
		DecModel*   model,   /**< model pointer */
		DspParams*  par,     /**< parameters */
		DspMessage* message, /**< message pointer */
		MPI_Comm    comm     /**< MPI communicator */):
DwSolverSerial(model, par, message), comm_(comm) {
	MPI_Comm_size(comm_, &comm_size_);
	MPI_Comm_rank(comm_, &comm_rank_);
}

DwSolverMpi::~DwSolverMpi() {
	comm_ = MPI_COMM_NULL;
}

DSP_RTN_CODE DwSolverMpi::init() {
	BGN_TRY_CATCH
	worker_ = new DwWorkerMpi(model_, par_, message_, comm_);
	if (comm_rank_ == 0) {
		master_ = new DwMasterTr(worker_);
		DSP_RTN_CHECK_THROW(master_->init());
	}
	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

DSP_RTN_CODE DwSolverMpi::solve() {
	BGN_TRY_CATCH
	if (comm_rank_ == 0) {
		DSP_RTN_CHECK_THROW(master_->solve());
	} else {
		DSP_RTN_CHECK_THROW(dynamic_cast<DwWorkerMpi*>(worker_)->receiver());
	}
	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

DSP_RTN_CODE DwSolverMpi::finalize() {
	BGN_TRY_CATCH
	if (comm_rank_ == 0)
		DSP_RTN_CHECK_THROW(master_->finalize());
	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

bool DwSolverMpi::chooseBranchingObjects(
		DspBranch*& branchingUp, /**< [out] branching-up object */
		DspBranch*& branchingDn  /**< [out] branching-down object */) {
	int branched = 0;
	if (comm_rank_ == 0)
		branched = master_->chooseBranchingObjects(branchingUp, branchingDn);
	return branched;
}

void DwSolverMpi::setBranchingObjects(const DspBranch* branchobj) {
	BGN_TRY_CATCH

	if (comm_rank_ == 0)
		master_->setBranchingObjects(branchobj);

	END_TRY_CATCH(;)
}
