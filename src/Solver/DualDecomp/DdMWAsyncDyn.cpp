/*
 * DdMWAsyncDyn.cpp
 *
 *  Created on: Feb 20, 2018
 *      Author: Kibaek Kim
 */

//#define DSP_DEBUG
//#define DSP_DEBUG_QUEUE

#define SOLUTION_KEY_TO_STOP -999

#include "Solver/DualDecomp/DdMasterAtr.h"
#include "Solver/DualDecomp/DdMWAsyncDyn.h"

DdMWAsyncDyn::DdMWAsyncDyn(
		MPI_Comm     comm,   /**< MPI communicator */
		DecModel *   model,  /**< model pointer */
		DspParams *  par,    /**< parameters */
		DspMessage * message /**< message pointer */):
DdMWAsync(comm, model, par, message) {
	/** nothignt to do */
}

DdMWAsyncDyn::DdMWAsyncDyn(const DdMWAsyncDyn& rhs) :
DdMWAsync(rhs) {}

DdMWAsyncDyn::~DdMWAsyncDyn() {
	/** nothignt to do */
}

bool DdMWAsyncDyn::chooseQueueElement(int& qid, double*& qsol, int& nsubprobs, int*& subindex) {
	/** initialize queue ID */
	qid = -1;

	for (unsigned k = 0; k < q_indicator_.size(); ++k) {
		int kk = k;
		/** LIFO? */
		if (par_->getBoolParam("DD/ASYNC/FIFO") == false)
			kk = q_indicator_.size() - 1 - k;
		for (int s = 0; s < model_->getNumSubproblems(); ++s) {
			if (q_indicator_[kk][s] == Q_NOT_ASSIGNED) {
				qid = q_id_[kk];
				qsol = q_solution_[kk];
				nsubprobs = 1;
				subindex[0] = s;
				q_indicator_[kk][s] = Q_ASSIGNED;
				break;
			}
		}
		if (qid > -1) break;
	}

	return (qid > -1);
}

DSP_RTN_CODE DdMWAsyncDyn::setWorkerLb(DdWorkerLB* workerlb, int nsubprobs, int* subindex, double* buf, double bestprimobj) {
	BGN_TRY_CATCH

	if (nsubprobs == 1) {
		workerlb->createProblem(1, subindex);
		workerlb->subprobs_[0]->theta_ = buf[0];
		workerlb->subprobs_[0]->updateProblem(buf + 1, bestprimobj);

		/** TODO: need to check if cut generation still works */
		/** apply Benders cuts */
		if (cutsToAdd_->sizeCuts() > 0 && (parFeasCuts_ >= 0 || parOptCuts_ >= 0)) {
			workerlb->subprobs_[0]->pushCuts(cutsToAdd_);
			DSPdebugMessage("Rank %d pushed %d Benders cuts.\n", comm_rank_, cutsToAdd_->sizeCuts());
		}

	} else if (nsubprobs > 1) 
		throw "Dynamic allocation: the worker received more than one subproblem.";

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

