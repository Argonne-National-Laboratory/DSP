/*
 * DwWorkerMpi.cpp
 *
 *  Created on: Nov 18, 2016
 *      Author: kibaekkim
 */

#include <DantzigWolfe/DwWorkerMpi.h>

DwWorkerMpi::DwWorkerMpi(
		DecModel * model,
		DspParams * par,
		DspMessage * message,
		MPI_Comm comm):
DwWorker(model, par, message),
comm_(comm) {
	MPI_Comm_rank(comm_, &comm_rank_);
	MPI_Comm_size(comm_, &comm_size_);

	/** calculate the size of piA vector */
	if (comm_rank_ > 0) {
		if (model_->isStochastic()) {
			TssModel* tss = dynamic_cast<TssModel*>(model_);
			npiA_ = tss->getNumScenarios() * (tss->getNumCols(0) + tss->getNumCols(1));
		} else
			npiA_ = model_->getNumSubproblemCouplingCols(0);
	}

	/** broadcast the size of piA vector */
	MPI_Bcast(&npiA_, 1, MPI_INT, 0, comm_);

	/** create a local vector memory */
	piA_ = new double [npiA_];

	/** get the number of subproblems */
	MPI_Reduce(&parProcIdxSize_, &nsubprobs_, 1, MPI_INT, MPI_SUM, 0, comm_);

	/** start coordinator */
	if (comm_rank_ > 0)
		coordinator();
}

DwWorkerMpi::~DwWorkerMpi() {
	FREE_ARRAY_PTR(piA_);
}

DSP_RTN_CODE DwWorkerMpi::coordinator() {
	int signal;
	std::vector<int> indices;
	std::vector<int> statuses;
	std::vector<double> cxs;
	std::vector<double> objs;
	std::vector<CoinPackedVector*> sols;
	bool terminate = false;
	CoinError::printErrors_ = true;

	BGN_TRY_CATCH

	while (terminate == false) {
		/** receive a signal */
		MPI_Bcast(&signal, 1, MPI_INT, 0, comm_);

		switch(signal) {
		case sig_generateCols:
			DSP_RTN_CHECK_RTN_CODE(
					generateCols(-1, NULL, indices, statuses, cxs, objs, sols));
			break;
		case sig_setColBounds:
			setColBounds(-1, 0.0, 0.0);
			break;
		case sig_terminate:
			terminate = true;
			break;
		default:
			CoinError("Unexpected signal", "DwWorkerMpi", "coordinator");
			break;
		}
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DwWorkerMpi::generateCols(
		int phase,                           /**< [in] phase of the master */
		const double* piA,                   /**< [in] piA */
		std::vector<int>& indices,           /**< [out] subproblem indices */
		std::vector<int>& statuses,          /**< [out] solution status */
		std::vector<double>& cxs,            /**< [out] solution times original objective coefficients */
		std::vector<double>& objs,           /**< [out] subproblem objective values */
		std::vector<CoinPackedVector*>& sols /**< [out] subproblem coupling column solutions */) {
#define FREE_MEMORY \
		FREE_ARRAY_PTR(_indices) \
		FREE_ARRAY_PTR(_statuses) \
		FREE_ARRAY_PTR(_cxs) \
		FREE_ARRAY_PTR(_objs)

	int* recvcounts = NULL;
	int* displs = NULL;
	int* _indices = NULL;
	int* _statuses = NULL;
	double* _cxs = NULL;
	double* _objs = NULL;
	std::vector<CoinPackedVector*> _sols;

	BGN_TRY_CATCH

	if (comm_rank_ == 0) {
		/** send signal */
		MPI_Bcast(&sig_generateCols, 1, MPI_INT, 0, comm_);

		/** copy piA vector */
		CoinCopyN(piA, npiA_, piA_);

		recvcounts = new int [comm_size_];
		displs = new int [comm_size_];
		_indices = new int [nsubprobs_];
		_statuses = new int [nsubprobs_];
		_cxs = new double [nsubprobs_];
		_objs = new double [nsubprobs_];
	}

	/** The root rank gathers the number of subproblems for each process. */
	MPI_Gather(&parProcIdxSize_, 1, MPI_INT, recvcounts, 1, MPI_INT, 0, comm_);
	/** calculate displacement of the receive buffer */
	displs[0] = 0;
	for (int i = 1; i < comm_size_; ++i)
		displs[i] = displs[i-1] + recvcounts[i-1];

	/** synchronize information */
	MPI_Bcast(&phase, 1, MPI_INT, 0, comm_);
	MPI_Bcast(piA_, npiA_, MPI_DOUBLE, 0, comm_);

	DSP_RTN_CHECK_RTN_CODE(
			DwWorker::generateCols(phase, piA, indices, statuses, cxs, objs, sols));

	/** synchronize information */
	MPI_Gatherv(&indices[0], indices.size(), MPI_INT,
			_indices, recvcounts, displs, MPI_INT, 0, comm_);
	MPI_Gatherv(&statuses[0], statuses.size(), MPI_INT,
			_statuses, recvcounts, displs, MPI_INT, 0, comm_);
	MPI_Gatherv(&cxs[0], cxs.size(), MPI_DOUBLE,
			_cxs, recvcounts, displs, MPI_DOUBLE, 0, comm_);
	MPI_Gatherv(&objs[0], objs.size(), MPI_DOUBLE,
			_objs, recvcounts, displs, MPI_DOUBLE, 0, comm_);
	MPIgatherCoinPackedVectors(comm_, _sols, sols);

	/** construct output arguments */
	if (comm_rank_ == 0) {
		indices.clear();
		statuses.clear();
		cxs.clear();
		objs.clear();
		sols.clear();
		for (int i = 0; i < nsubprobs_; ++i) {
			indices.push_back(_indices[i]);
			statuses.push_back(_statuses[i]);
			cxs.push_back(_cxs[i]);
			objs.push_back(_objs[i]);
			sols.push_back(_sols[i]);
			_sols[i] = NULL;
		}
		_sols.clear();
	}

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
}

void DwWorkerMpi::setColBounds(int j, double lb, double ub) {

	BGN_TRY_CATCH

	/** send signal */
	if (comm_rank_ == 0)
		MPI_Bcast(&sig_setColBounds, 1, MPI_INT, 0, comm_);

	/** synchronize information */
	MPI_Bcast(&j, 1, MPI_INT, 0, comm_);
	MPI_Bcast(&lb, 1, MPI_DOUBLE, 0, comm_);
	MPI_Bcast(&ub, 1, MPI_DOUBLE, 0, comm_);

	DwWorker::setColBounds(j, lb, ub);

	END_TRY_CATCH(;)
}
