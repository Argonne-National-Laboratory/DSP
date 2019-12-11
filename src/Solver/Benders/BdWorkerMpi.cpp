/*
 * BdWorkerMpi.cpp
 *
 *  Created on: Sep 19, 2016
 *      Author: kibaekkim
 */

#include "Benders/BdWorkerMpi.h"

BdWorkerMpi::BdWorkerMpi(DecModel * model, DspParams * par, DspMessage * message, MPI_Comm comm) :
		BdWorker(model, par, message),
		comm_(comm),
		nx_(0),
		naux_(0),
		x_(NULL),
		status_(NULL),
		cut_index_(NULL),
		recvcounts_(NULL),
		displs_(NULL) {
	MPI_Comm_rank(comm_, &comm_rank_);
	MPI_Comm_size(comm_, &comm_size_);
}

BdWorkerMpi::~BdWorkerMpi() {
	FREE_ARRAY_PTR(x_);
	FREE_ARRAY_PTR(status_);
	FREE_ARRAY_PTR(cut_index_);
	FREE_ARRAY_PTR(recvcounts_);
	FREE_ARRAY_PTR(displs_);
}

int BdWorkerMpi::init() {
	BGN_TRY_CATCH

	/** Receive number of columns */
	MPI_Bcast(&nx_, 1, MPI_INT, 0, comm_);
	/** Receive number of auxiliary columns */
	MPI_Bcast(&naux_, 1, MPI_INT, 0, comm_);
	/** allocate memory for the solution */
	x_ = new double [nx_];
	if (comm_ == 0) {
		status_    = new int [parProcIdxSize_];
		cut_index_ = new int [parProcIdxSize_];
		recvcounts_ = new int [comm_size_];
		displs_     = new int [comm_size_];
		/** Subproblems are assigned to each process in round-and-robin fashion. */
		for (int i = 0; i < comm_size_; ++i) {
			recvcounts_[i] = 0;
			for (int s = i; s < parProcIdxSize_; s += comm_size_) {
				recvcounts_[i]++;
				cut_index_[s] = i;
			}
			displs_[i] = i == 0 ? 0 : displs_[i-1] + recvcounts_[i-1];
		}
	}

	/** run non-root processors */
	if (comm_ != 0) {
		OsiCuts cs;
		generateCuts(nx_, naux_, NULL, cs);
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE BdWorkerMpi::generateCuts(int nx, int naux, const double* x, OsiCuts& cs) {

#define FREE_MEMORY \
		FREE_2D_ARRAY_PTR(cutval, parProcIdxSize_) \
		FREE_ARRAY_PTR(cutrhs)

	double** cutval = NULL;
	double*  cutrhs = NULL;
	int message = 1;

	BGN_TRY_CATCH

	while (comm_ != 0) {
		/** Receive message */
		MPI_Bcast(&message, 1, MPI_INT, 0, comm_);

		/** TODO: make sure the master sends 0 at the end of solution. */
		/** done if message == 0; go, otherwise. */
		if (message == 0) break;

		/** Receive master solution */
		MPI_Bcast(x_, nx_, MPI_DOUBLE, 0, comm_);

		/** allocate memory for cuts */
		cutval = new double * [parProcIdxSize_];
		cutrhs = new double [parProcIdxSize_];
		for (int s = 0; s < parProcIdxSize_; ++s)
			cutval[s] = NULL;
		CoinZeroN(cutrhs, parProcIdxSize_);

		/** generate cuts */
		bdsub_->generateCuts(nx_, x_, cutval, cutrhs);

		/** collect cuts */
		collectCuts(nx_, naux_, cutval, cutrhs, cs);
	}

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

DSP_RTN_CODE BdWorkerMpi::collectCuts(
		int nx,
		int naux,
		double** cut,
		double* rhs, OsiCuts& cs) {
#define FREE_MEMORY \
	FREE_2D_ARRAY_PTR(aggcut, naux) \
	FREE_ARRAY_PTR(aggrhs)

	double** aggcut = NULL;
	double* aggrhs = NULL;
	CoinPackedVector vec;
	OsiCuts cuts;

	BGN_TRY_CATCH

	/** construct cut information for each process */
	for (int i = 0; i < bdsub_->getNumSubprobs(); ++i) {
		/** initialize vector */
		vec.clear();
		/** set it as sparse */
		for (int j = 0; j < nx; ++j) {
			if (fabs(cut[i][j]) > 1e-10)
				vec.insert(j, cut[i][j]);
		}

		/** create row cut */
		OsiRowCut rc;
		rc.setRow(vec);
		rc.setUb(COIN_DBL_MAX);
		rc.setLb(rhs[i]);

		/** add cut */
		cuts.insert(rc);
	}

	/** Send cut generation status to the master */
	MPI_Gatherv(bdsub_->getStatuses(), bdsub_->getNumSubprobs(), MPI_INT,
			status_, recvcounts_, displs_, MPI_INT, 0, comm_);

	/** Send cuts to the master */
	OsiCuts tempcuts;
	MPIgatherOsiCuts(comm_, cuts, tempcuts);

	/** cleanup cuts */
	for (int i = 0; i < cuts.sizeCuts(); ++i) {
		OsiRowCut * rc = cuts.rowCutPtr(i);
		FREE_PTR(rc);
	}
	cuts.dumpCuts();

	/** the root generates cuts; the others are done. */
	if (comm_ == 0) {
		/** is there a feasibility cut? */
		int fcut = -1;
		for (int i = 0; i < model_->getNumSubproblems(); ++i)
			if (status_[i] == DSP_STAT_PRIM_INFEASIBLE) {
				fcut = i;
				break;
			}

		if (fcut > -1) {
			/** create row cut */
			OsiRowCut rc = tempcuts.rowCut(fcut);
			/** store cut */
			cs.insert(rc);
		} else {
			/** total number of subproblems */
			int nsubprobs = model_->getNumSubproblems();

			/** allocate memory for cuts */
			aggcut = new double * [naux];
			aggrhs = new double [naux];
			for (int s = 0; s < naux; ++s) {
				aggcut[s] = new double [nx];
				CoinZeroN(aggcut, nx);
			}
			CoinZeroN(aggrhs, naux);

			/** aggregate cuts */
			for (int s = 0; s < nsubprobs; ++s) {
				OsiRowCut* rc = tempcuts.rowCutPtr(s);
				const CoinPackedVector row = rc->row();
				int i = cut_index_[s] % naux;
				for (int j = 0; j < row.getNumElements(); ++j)
					aggcut[i][row.getIndices()[j]] += row.getElements()[j];
				aggrhs[i] += rc->lb();
			}

			/** construct cuts */
			for (int i = 0; i < naux; ++i) {

				/** initialize vector */
				vec.clear();

				/** set it as sparse */
				for (int j = 0; j < nx; ++j) {
					if (fabs(aggcut[i][j]) > 1e-10)
						vec.insert(j, aggcut[i][j]);
				}
				vec.insert(nx + i, 1.0);

				/** create row cut */
				OsiRowCut rc;
				rc.setRow(vec);
				rc.setUb(COIN_DBL_MAX);
				rc.setLb(aggrhs[i]);

				/** store cut */
				cs.insert(rc);
			}
		}
	}

	for (int i = 0; i < tempcuts.sizeCuts(); ++i) {
		OsiRowCut * rc = tempcuts.rowCutPtr(i);
		FREE_PTR(rc);
	}
	tempcuts.dumpCuts();

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}
