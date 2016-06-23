/*
 * DdMWSerial.cpp
 *
 *  Created on: Apr 8, 2016
 *      Author: kibaekkim
 */

#include <Solver/DualDecomp/DdMWSerial.h>

DdMWSerial::DdMWSerial(): DdMW(){
	// TODO Auto-generated constructor stub
}

DdMWSerial::DdMWSerial(
		DdMaster *        master, /**< master problem */
		vector<DdWorker*> worker  /**< worker for finding lower bounds */):
DdMW(master, worker){
	// TODO Auto-generated constructor stub

}

DdMWSerial::~DdMWSerial() {
	// TODO Auto-generated destructor stub
}

DSP_RTN_CODE DdMWSerial::run() {

	BGN_TRY_CATCH

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMWSerial::init() {
}

DSP_RTN_CODE DdMWSerial::runMaster() {

	BGN_TRY_CATCH

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMWSerial::runWorker() {

	DdMasterSync * master = NULL;
	DdWorkerLB * workerlb = NULL;
	DdWorkerUB * workerub = NULL;

	BGN_TRY_CATCH

	/** retrieve DdMasterSync */
	master = dynamic_cast<DdMasterSync*>(master_);

	/** retrieve DdWorkerLB and DdWorkerUB */
	assert(worker_.size()>0);
	assert(worker_[0]->getType()==DdWorker::LB);
	workerlb = dynamic_cast<DdWorkerLB*>(worker_[0]);
	if (worker_.size() > 1)
	{
		assert(worker_[1]->getType()==DdWorker::UB);
		workerub = dynamic_cast<DdWorkerUB*>(worker_[1]);
	}

	double * master_primsol = const_cast<double*>(master_->getPrimalSolution());
	for (unsigned s = 0; s < workerlb->subprobs_.size(); ++s) {
		workerlb->subprobs_[s]->theta_ = master_primsol[j];
		workerlb->subprobs_[s]->updateProblem(master_primsol + model_->getNumSubproblems());
	}

	/** Solve subproblems assigned to each process  */
	workerlb->solve();

	/** create send buffer */
	for (unsigned s = 0, pos = 0; s < workerlb->subprobs_.size(); ++s)
	{
		sendbuf[pos++] = static_cast<double>(workerlb->subprobs_[s]->sind_);
		sendbuf[pos++] = workerlb->subprobs_[s]->getPrimalBound();
		sendbuf[pos++] = workerlb->subprobs_[s]->getDualBound();
		CoinCopyN(workerlb->subprobs_[s]->si_->getSolution(), workerlb->subprobs_[s]->ncols_coupling_, sendbuf + pos);
		pos += model_->getNumSubproblemCouplingCols(workerlb->subprobs_[s]->sind_);
		message_->print(5, "MW -> worker %d, subprob %d primobj %+e dualobj %+e\n",
				comm_rank_, workerlb->subprobs_[s]->sind_,
				workerlb->subprobs_[s]->getPrimalBound(), workerlb->subprobs_[s]->getDualBound());
	}

	/** calculate upper bounds */
	if (workerub & solutions.size() > 0) {
		upperbounds.clear();
		for (unsigned i = 0; i < solutions.size(); ++i) {
			/** fix coupling solutions */
			workerub->fixCouplingVariableValues(solutions[i]);

			/** solve */
			workerub->solve();

			/** take minimum objective */
			double sumprimobj = 0.0;
			for (unsigned s = 0; s < workerub->subprobs_.size(); ++s)
				sumprimobj += workerub->subprobs_[s]->getPrimalBound();
			upperbounds.push_back(sumprimobj);
		}
	}

	for (unsigned i = 0; i < solutions.size(); ++i)
		FREE_PTR(solutions[i]);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMWSerial::finalize() {
}
