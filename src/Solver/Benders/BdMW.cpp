/*
 * BdMW.cpp
 *
 *  Created on: Feb 17, 2016
 *      Author: kibaekkim
 */

#include <algorithm>
#include "Solver/Benders/BdMW.h"

BdMW::BdMW(MPI_Comm comm, BdMaster * master, BdWorker * worker):
	BaseMasterWorker(comm),
	master_(master), worker_(worker), primsol_(NULL) {}

BdMW::~BdMW()
{
	FREE_ARRAY_PTR(primsol_);
}

STO_RTN_CODE BdMW::run()
{
	BGN_TRY_CATCH

	/** initialize */
	init();

	/** run master process */
	runMaster();

	/** run worker processes */
	runWorker();

	/** finalize */
	finalize();

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
}

STO_RTN_CODE BdMW::init()
{
	BGN_TRY_CATCH

	if (comm_rank_ == 0)
		primsol_  = new double [master_->getModelPtr()->getFullModelNumCols()];

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
}

STO_RTN_CODE BdMW::runMaster()
{
	if (!master_)
		return STO_RTN_OK;

	BGN_TRY_CATCH

	/** solve */
	master_->solve();

	/** Tell workers we are done */
	int message = MASTER_STOPPED;
	MPI_Bcast(&message, 1, MPI_INT, 0, comm_);

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
}

STO_RTN_CODE BdMW::runWorker()
{
#define FREE_MEMORY \
	FREE_ARRAY_PTR(solution); \
	FREE_2D_ARRAY_PTR(model->getNumSubproblems(),cutval); \
	FREE_ARRAY_PTR(status); \
	FREE_ARRAY_PTR(cutrhs);

	if (!worker_)
		return STO_RTN_OK;

	OsiCuts cuts, tempcuts;
	int message;
	int ncols;
	double * solution = NULL;
	double ** cutval = NULL;
	double * cutrhs = NULL;
	CoinPackedVector vec;
	int * status = NULL;
	DecModel * model = NULL;
	DspParams * par = NULL;
	BdSub * bdsub = NULL;

	BGN_TRY_CATCH

	/** internal pointers */
	model = worker_->getModelPtr();
	par   = worker_->getParPtr();
	bdsub = worker_->getBdSubPtr();

	int parProcIdxSize = par->getIntPtrParamSize("ARR_PROC_IDX");

	ncols = model->getNumSubproblemCouplingCols(0) + par->getIntParam("BD/NUM_CUTS_PER_ITER");
	solution = new double [ncols];
	cutval = new double * [model->getNumSubproblems()];
	cutrhs = new double [model->getNumSubproblems()];
	status = new int [parProcIdxSize];

	/** Wait for message from the master */
	MPI_Bcast(&message, 1, MPI_INT, 0, comm_);
	DSPdebugMessage("[%d]: Received message [%d]\n", comm_rank_, message);

	/** Parse the message */
	while (message == MASTER_NEEDS_CUTS)
	{
		/** Receive master solution */
		MPI_Bcast(solution, ncols, MPI_DOUBLE, 0, comm_);

		/** Generate cuts */
		bdsub->generateCuts(ncols, solution, cutval, cutrhs);
		for (int s = 0, ss = 0; s < bdsub->getNumSubprobs(); ++s)
		{
			/** initialize vector */
			vec.clear();

			/** set it as sparse */
			for (int j = 0; j < ncols; ++j)
				if (fabs(cutval[s][j]) > 1e-10)
					vec.insert(j, cutval[s][j]);

			/** free memory */
			FREE_ARRAY_PTR(cutval[s]);

			if (fabs(cutrhs[s]) < 1e-10)
				cutrhs[s] = 0.0;

			OsiRowCut rc;
			rc.setRow(vec);
			rc.setUb(COIN_DBL_MAX); /** TODO: for minimization */
			rc.setLb(cutrhs[s]);

			//DSPdebug(rc.print());
			cuts.insert(rc);

			/** get status */
			status[ss++] = bdsub->getStatus(s);
			DSPdebugMessage("[%d]: status[%d] %d\n", comm_rank_, ss-1, status[ss-1]);
		}
		DSPdebugMessage("[%d]: Found %d cuts\n", comm_rank_, cuts.sizeCuts());

		/** Send cut generation status to the master */
		MPI_Gatherv(status, parProcIdxSize, MPI_INT, NULL, NULL, NULL, MPI_INT, 0, comm_);

		/** Send cuts to the master */
		MPIgatherOsiCuts(comm_, cuts, tempcuts);

		/** cleanup cuts */
		for (int i = 0; i < cuts.sizeCuts(); ++i)
		{
			OsiRowCut * rc = cuts.rowCutPtr(i);
			FREE_PTR(rc);
		}
		cuts.dumpCuts();

		/** Wait for message from the master */
		MPI_Bcast(&message, 1, MPI_INT, 0, comm_);
		DSPdebugMessage("[%d]: Received message [%d]\n", comm_rank_, message);
	}

	END_TRY_CATCH_RTN(FREE_MEMORY,STO_RTN_ERR)

	FREE_MEMORY

	return STO_RTN_OK;
#undef FREE_MEMORY
}

STO_RTN_CODE BdMW::finalize()
{
#define FREE_MEMORY            \
	FREE_ARRAY_PTR(nsubprobs)  \
	FREE_ARRAY_PTR(subindices) \
	FREE_ARRAY_PTR(displs)     \
	FREE_ARRAY_PTR(sizesols)   \
	FREE_ARRAY_PTR(sizesolsp)  \
	FREE_ARRAY_PTR(subsols)

	int * nsubprobs  = NULL;
	int * subindices = NULL;
	int * displs     = NULL;
	int * sizesols   = NULL; /**< size of solution vectors for each subproblem */
	int * sizesolsp  = NULL; /**< size of solution vectors for each process */
	double * subsols = NULL;

	BGN_TRY_CATCH

	if (comm_rank_ == 0)
	{
		nsubprobs  = new int [comm_size_];
		subindices = new int [master_->getModelPtr()->getNumSubproblems()];
		displs     = new int [comm_size_];
		sizesols   = new int [master_->getModelPtr()->getNumSubproblems()];
		sizesolsp  = new int [comm_size_];

		/** msg#1. receive number of subproblems */
		MPI_Gather(NULL, 0, MPI_INT, nsubprobs, comm_size_, MPI_INT, 0, comm_);

		/** set displacement */
		displs[0] = 0;
		for (int i = 1; i < comm_size_; ++i)
			displs[i] = displs[i-1] + nsubprobs[i-1];

		/** msg#2. receive subproblem indices */
		MPI_Gatherv(NULL, 0, MPI_INT, subindices, nsubprobs, displs, MPI_INT, 0, comm_);

		/** msg#3. receive size of solution vectors for each subproblem */
		MPI_Gatherv(NULL, 0, MPI_INT, sizesols, nsubprobs, displs, MPI_INT, 0, comm_);

		/** msg#4. receive size of solution vectors for each process */
		MPI_Gather(NULL, 0, MPI_INT, sizesolsp, comm_size_, MPI_INT, 0, comm_);

		/** set displacement */
		displs[0] = 0;
		for (int i = 1; i < comm_size_; ++i)
			displs[i] = displs[i-1] + sizesolsp[i-1];
		subsols = new double [displs[comm_size_-1]+sizesolsp[comm_size_-1]];

		/** msg#5. receive subproblem solutions */
		MPI_Gatherv(NULL, 0, MPI_DOUBLE, subsols, sizesolsp, displs, MPI_DOUBLE, 0, comm_);

		/** copy master solution */
		CoinCopyN(master_->getPrimalSolution(), master_->getModelPtr()->getNumSubproblemCouplingCols(0), primsol_);

		vector<double*> vsubsols(master_->getModelPtr()->getNumSubproblems());
		vector<int> vsubsollens(master_->getModelPtr()->getNumSubproblems());
		for (int i = 0, j = 0; i < master_->getModelPtr()->getNumSubproblems(); ++i)
		{
			vsubsols[subindices[i]] = subsols + j;
			vsubsollens[subindices[i]] = sizesols[i];
			j += sizesols[i];
		}

		/** copy subproblem solution */
		int pos = master_->getModelPtr()->getNumSubproblemCouplingCols(0);
		for (int i = 0; i < master_->getModelPtr()->getNumSubproblems(); ++i)
		{
			CoinCopyN(vsubsols[i], vsubsollens[i], primsol_ + pos);
			pos += vsubsollens[i];
		}
	}
	else
	{
		int num = worker_->getBdSubPtr()->getNumSubprobs();
		sizesols = new int [num];

		/** msg#1. send number of subproblems */
		MPI_Gather(&num, 1, MPI_INT, NULL, 0, MPI_INT, 0, comm_);

		/** mgs#2. send subproblem indices */
		MPI_Gatherv(worker_->getBdSubPtr()->getSubprobIndices(), worker_->getBdSubPtr()->getNumSubprobs(),
				MPI_INT, NULL, NULL, NULL, MPI_INT, 0, comm_);

		int scount = 0;
		for (int i = 0; i < worker_->getBdSubPtr()->getNumSubprobs(); ++i)
		{
			sizesols[i] = worker_->getBdSubPtr()->getNumCols(i);
			scount += worker_->getBdSubPtr()->getNumCols(i);
		}

		/** msg#3. send size of solution vectors for each subproblem*/
		MPI_Gatherv(sizesols, worker_->getBdSubPtr()->getNumSubprobs(), MPI_INT, NULL, NULL, NULL, MPI_INT, 0, comm_);

		/** msg#4. send size of solution vectors for each process */
		MPI_Gather(&scount, 1, MPI_INT, NULL, 0, MPI_INT, 0, comm_);

		subsols = new double [scount];
		for (int i = 0, pos = 0; i < worker_->getBdSubPtr()->getNumSubprobs(); ++i)
		{
			CoinCopyN(worker_->getBdSubPtr()->getSolution(i),
					worker_->getBdSubPtr()->getNumCols(i),
					subsols + pos);
			pos += worker_->getBdSubPtr()->getNumCols(i);
		}

		/** msg#5. send subproblem solutions */
		MPI_Gatherv(subsols, scount, MPI_DOUBLE, NULL, NULL, NULL, MPI_DOUBLE, 0, comm_);
	}

	END_TRY_CATCH_RTN(FREE_MEMORY,STO_RTN_ERR)

	FREE_MEMORY

	return STO_RTN_OK;
#undef FREE_MEMORY
}
