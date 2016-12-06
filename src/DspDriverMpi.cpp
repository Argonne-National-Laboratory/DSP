/*
 * DspDriverMpi.cpp
 *
 *  Created on: Dec 5, 2016
 *      Author: kibaekkim
 */

//#define DSP_DEBUG
#include <DspDriverMpi.h>
#include <DantzigWolfe/DwSolverMpi.h>

DspDriverMpi::DspDriverMpi(
		DecModel* model, /**< model pointer */
		DspParams* par,  /**< parameter pointer */
		MPI_Comm comm    /**< MPI Communicator */):
DspDriver(model, par), comm_(comm) {
	MPI_Comm_size(comm_, &comm_size_);
	MPI_Comm_rank(comm_, &comm_rank_);
	DSPdebugMessage("Rank %d started.\n", comm_rank_);
}

DSP_RTN_CODE DspDriverMpi::init() {
	BGN_TRY_CATCH

	solver_ = new DwSolverMpi(model_, par_, message_, comm_);
	DSP_RTN_CHECK_RTN_CODE(solver_->init());

	if (comm_rank_ == 0) {
		/** create an Alps model */
		alps_ = new DspModel(solver_);
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DspDriverMpi::run() {
	BGN_TRY_CATCH
	if (comm_rank_ == 0) {
		DSP_RTN_CHECK_RTN_CODE(alps_->solve());
		/** send signal */
		int sig = DwWorkerMpi::sig_terminate;
		MPI_Bcast(&sig, 1, MPI_INT, 0, comm_);
		DSPdebugMessage("Rank 0 sent signal %d.\n", sig);
	} else {
		DSP_RTN_CHECK_RTN_CODE(solver_->solve());
	}
	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}
