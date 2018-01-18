/*
 * DwWorkerPips.cpp
 *
 *  Created on: Dec 15, 2016
 *      Author: Kibaek Kim
 */

//#define DSP_DEBUG
#include "Model/TssModel.h"
#include "Solver/DantzigWolfe/DwWorkerPips.h"
#include "Utility/DspMpi.h"

DwWorkerPips::DwWorkerPips(
		DecModel * model,
		DspParams * par,
		DspMessage * message,
		MPI_Comm comm):
DwWorkerMpi(model, par, message, comm), pips_(NULL) {
}

DSP_RTN_CODE DwWorkerPips::receiver() {
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
		DSPdebugMessage("Rank %d received signal %d.\n", comm_rank_, signal);

		switch(signal) {
		case sig_generateCols:
			DSP_RTN_CHECK_RTN_CODE(
					generateCols(-1, NULL, indices, statuses, cxs, objs, sols));
			indices.clear();
			statuses.clear();
			cxs.clear();
			objs.clear();
			sols.clear();
			break;
		case sig_initPips:
			initPips();
			break;
		case sig_solvePips:
			DSP_RTN_CHECK_RTN_CODE(solvePips());
			break;
		case sig_clearMats:
			clearMats();
			break;
		case sig_setColBounds:
			setColBounds(0, NULL, NULL, NULL);
			break;
		case sig_terminate:
			terminate = true;
			break;
		default:
			CoinError("Unexpected signal", "DwWorkerPips", "coordinator");
			break;
		}
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

void DwWorkerPips::initPips(int nscen, int ncols) {

	if (comm_rank_ == 0) {
		/** send signal */
		int sig = sig_initPips;
		DSPdebugMessage("Rank %d sends signal %d.\n", comm_rank_, sig);
		MPI_Bcast(&sig, 1, MPI_INT, 0, comm_);
	}
	MPI_Bcast(&nscen, 1, MPI_INT, 0, comm_);
	MPI_Bcast(&ncols, 1, MPI_INT, 0, comm_);

	//printf("Rank %d creates PipsInterface(%d, %d).\n", comm_rank_, nscen, ncols);
	pips_ = new PipsInterface(nscen, ncols);

	if (comm_rank_ == 0) {
		/** initialize rstart for adding rows to workers */
		rstart_.resize(nscen, ncols);
	}
}

DSP_RTN_CODE DwWorkerPips::solvePips(double weight) {
#define FREE_MEMORY \
	for (int i = 0; i < vecs.size(); ++i) \
		FREE_PTR(vecs[i]);

	std::vector<CoinPackedVector*> vecs;
	std::vector<double> rlbd, rubd;
	int nrows;
	std::vector<int> nscen(comm_size_), scens, rcnt(comm_size_), displs(comm_size_);

	BGN_TRY_CATCH

	if (comm_rank_ == 0) {
		/** send signal */
		int sig = sig_solvePips;
		DSPdebugMessage("Rank %d sends signal %d.\n", comm_rank_, sig);
		MPI_Bcast(&sig, 1, MPI_INT, 0, comm_);
	}

	/** time lapse for communication */
	double time_to_comm = MPI_Wtime();

	/** gather number of subproblems for each rank */
	MPI_Gather(&parProcIdxSize_, 1, MPI_INT, &nscen[0], 1, MPI_INT, 0, comm_);
	if (comm_rank_ == 0) {
		for (int j = 0; j < comm_size_; ++j) {
			rcnt[j] = nscen[j];
			if (j == 0)
				displs[j] = 0;
			else
				displs[j] = displs[j-1] + rcnt[j-1];
		}
		scens.resize(displs[comm_size_-1] + rcnt[comm_size_-1]);
	}
	MPI_Gatherv(parProcIdx_, parProcIdxSize_, MPI_INT, &scens[0], &rcnt[0], &displs[0], MPI_INT, 0, comm_);

	/** reset (clear matrices in) PIPS object */
	//pips_->clearMatricesTW();

	if (comm_rank_ == 0) {
#if 0
		//printf("Distributing data...\n");
		for (int j = 0; j < comm_size_; ++j) {
			//printf("comm_rank %d: nscen %d: ", j, nscen[j]);
			for (int i = 0; i < nscen[j]; ++i)
				printf("%d ", scens[displs[j]+i]);
			printf("\n");
		}
#endif
		pips_->retrieveBundles(rstart_, vecs, rlbd, rubd, scens);
		nrows = rlbd.size();
		MPIbcastCoinPackedVectors(comm_, vecs);
		MPI_Bcast(&nrows, 1, MPI_INT, 0, comm_);
		MPI_Bcast(&rlbd[0], nrows, MPI_DOUBLE, 0, comm_);
		MPI_Bcast(&rubd[0], nrows, MPI_DOUBLE, 0, comm_);
		MPI_Bcast(&scens[0], nrows, MPI_INT, 0, comm_);
		
#if 0
		for (int i = 0; i < nrows; ++i) {
			printf("rlbd %e rubd %e, vec[%d]:\n", rlbd[i], rubd[i], i);
			DspMessage::printArray(vecs[i]);
		}
#endif
	} else {
		vecs.clear();
		MPIbcastCoinPackedVectors(comm_, vecs);
		MPI_Bcast(&nrows, 1, MPI_INT, 0, comm_);
		rlbd.resize(nrows); rubd.resize(nrows); scens.resize(nrows);
		MPI_Bcast(&rlbd[0], nrows, MPI_DOUBLE, 0, comm_);
		MPI_Bcast(&rubd[0], nrows, MPI_DOUBLE, 0, comm_);
		MPI_Bcast(&scens[0], nrows, MPI_INT, 0, comm_);
		for (int i = 0; i < nrows; ++i) {
			/** adjust indices of vector */
			vecs[i]->getIndices()[0] = scens[i];
			for (int j = 1; j < vecs[i]->getNumElements(); ++j)
				vecs[i]->getIndices()[j] += pips_->getNumScenarios() - 1;
#if 0
			if (comm_rank_ == 1) {
				printf("row vectori[%d] scens %d\n", i, scens[i]);
				DspMessage::printArray(vecs[i]);
			}
#endif
			pips_->addRow(*(vecs[i]), rlbd[i], rubd[i]);
		}
	}

	//printf("Broadcasting weight...\n");
	MPI_Bcast(&weight, 1, MPI_DOUBLE, 0, comm_);

#if 0
	for (int i = 0; i < comm_size_; ++i) {
		if (i == comm_rank_) {
			printf("comm rank %d\n", comm_rank_);
			pips_->print();
		}
		MPI_Barrier(comm_);
	}
#endif
	//if (comm_rank_ == 0) pips_->print();
	if (comm_rank_ == 0) printf("Time to distribute PIPS data: %.2f seconds.\n", MPI_Wtime() - time_to_comm);

	/** run PIPS */
	pips_->solve(weight);

	/** gather number of scenarios distributed in PIPS */
	int nscen_in_pips = pips_->scenarios_.size();
	std::vector<int> nscen_distributed_in_pips(comm_size_);
	MPI_Gather(&nscen_in_pips, 1, MPI_INT, &nscen_distributed_in_pips[0], 1, MPI_INT, 0, comm_);

	/** gather second-stage solution (theta) */
	std::vector<int> scen_distributed_in_pips(model_->getNumSubproblems());
	std::vector<double> thetas_distributed_in_pips(model_->getNumSubproblems());
	displs.resize(comm_size_, 0);
	for (int i = 1; i < comm_size_; ++i) 
		displs[i] = displs[i-1] + nscen_distributed_in_pips[i-1];
	MPI_Gatherv(&pips_->scenarios_[0], nscen_in_pips, MPI_INT, &scen_distributed_in_pips[0], &nscen_distributed_in_pips[0], &displs[0], MPI_INT, 0, comm_);
	MPI_Gatherv(&pips_->thetas_[0], nscen_in_pips, MPI_DOUBLE, &thetas_distributed_in_pips[0], &nscen_distributed_in_pips[0], &displs[0], MPI_DOUBLE, 0, comm_);
	if (comm_rank_ == 0)
		for (unsigned j = 0; j < scen_distributed_in_pips.size(); ++j)
			pips_->solution_[ scen_distributed_in_pips[ j ] ] = thetas_distributed_in_pips[ j ];
#if 0
	if (comm_rank_ == 0) {
		//printf("nscen_distributed_in_pips:\n");
		//DspMessage::printArray(nscen_distributed_in_pips.size(), &nscen_distributed_in_pips[0]);
		//printf("scen_distributed_in_pips:\n");
		//DspMessage::printArray(scen_distributed_in_pips.size(), &scen_distributed_in_pips[0]);
		printf("thetas_distributed_in_pips:\n");
		DspMessage::printArray(thetas_distributed_in_pips.size(), &thetas_distributed_in_pips[0]);
		printf("DwWorkerPips solution (%u):\n", pips_->solution_.size());
		DspMessage::printArray(pips_->solution_.size(), &pips_->solution_[0]);
	}
#endif

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

void DwWorkerPips::clearMats() {
	if (comm_rank_ == 0) {
		/** send signal */
		int sig = sig_clearMats;
		DSPdebugMessage("Rank %d sends signal %d.\n", comm_rank_, sig);
		MPI_Bcast(&sig, 1, MPI_INT, 0, comm_);
	}

	pips_->clearMatricesTW();

	if (comm_rank_ == 0)
		rstart_.resize(pips_->input_->nScenarios(), pips_->input_->nFirstStageVars());
}

