/*
 * DdMW.cpp
 *
 *  Created on: Feb 9, 2016
 *      Author: kibaekkim
 */

#include "Solver/DualDecomp/DdMW.h"

DdMW::DdMW(
		MPI_Comm comm,
		DdMaster * master,
		DdWorker * worker):
	BaseMasterWorker(comm),
	master_(master), worker_(worker),
	nsubprobs_(NULL), subprob_indices_(NULL), subprob_displs_(NULL),
	iteration_limit_(COIN_INT_MAX) {}

DdMW::~DdMW()
{
	FREE_ARRAY_PTR(nsubprobs_);
	FREE_ARRAY_PTR(subprob_indices_);
	FREE_ARRAY_PTR(subprob_displs_);
}

STO_RTN_CODE DdMW::run()
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

/** initialize */
STO_RTN_CODE DdMW::init()
{
	BGN_TRY_CATCH

	int   narrprocidx = 0;    /**< size of pointer parameter ARR_PROC_IDX */
	int * arrprocidx  = NULL; /**< pointer parameter ARR_PROC_IDX */

	if (comm_rank_ == 0)
	{
		narrprocidx     = master_->getParPtr()->getIntPtrParamSize("ARR_PROC_IDX");
		arrprocidx      = master_->getParPtr()->getIntPtrParam("ARR_PROC_IDX");
		nsubprobs_      = new int [comm_size_];
		subprob_displs_ = new int [comm_size_];
	}
	else
	{
		narrprocidx = worker_->getParPtr()->getIntPtrParamSize("ARR_PROC_IDX");
		arrprocidx  = worker_->getParPtr()->getIntPtrParam("ARR_PROC_IDX");
	}

	/** gather number of subproblems for each process */
	MPI_Gather(&narrprocidx, 1, MPI_INT, nsubprobs_, 1, MPI_INT, 0, comm_);

	if (comm_rank_ == 0)
	{
		int sum_of_nsubprobs = 0;
		for (int i = 0; i < comm_size_; ++i)
		{
			sum_of_nsubprobs += nsubprobs_[i];
			subprob_displs_[i] = i == 0 ? 0 : subprob_displs_[i-1] + nsubprobs_[i-1];
		}
		subprob_indices_ = new int [sum_of_nsubprobs];
	}

	/** gather subproblem indices for each process */
	MPI_Gatherv(arrprocidx, narrprocidx, MPI_INT, subprob_indices_, nsubprobs_, subprob_displs_, MPI_INT, 0, comm_);

	if (comm_rank_ == 0)
		sync_ = !master_->getParPtr()->getBoolParam("DD/ASYNC");
	else
		sync_ = !worker_->getParPtr()->getBoolParam("DD/ASYNC");

	/** free memory */
	arrprocidx = NULL;

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
}

/** run master process */
STO_RTN_CODE DdMW::runMaster()
{
	if (comm_rank_ != 0)
		return STO_RTN_OK;

	/** retrieve parameters */
	iteration_limit_ = master_->getParPtr()->getIntParam("ITER_LIM");

	if (sync_)
		return runMasterSync();
	else
	{
		/** need to run the first iteration of the synchronization,
		 * otherwise the master will become dual infeasible.
		 */

		master_->getMessagePtr()->print(1, "Running one iteration to collect bundles from all subproblems.\n");

		int temp_loglevel        = master_->getMessagePtr()->logLevel_;
		int temp_iteration_limit = iteration_limit_;
		master_->getMessagePtr()->logLevel_ = -1;
		iteration_limit_ = 0;

		runMasterSync();

		master_->getMessagePtr()->logLevel_ = temp_loglevel;
		iteration_limit_ = temp_iteration_limit;

		/** spawn the current solution to primsol_to_worker_ */
		master_->setInitSolution(master_->primsol_);

		return runMasterAsync();
	}
}

/** run master process */
STO_RTN_CODE DdMW::runMasterSync()
{
#define FREE_MEMORY \
	FREE_ARRAY_PTR(sendbuf) \
	FREE_ARRAY_PTR(scounts) \
	FREE_ARRAY_PTR(sdispls) \
	FREE_ARRAY_PTR(recvbuf) \
	FREE_ARRAY_PTR(rcounts) \
	FREE_ARRAY_PTR(rdispls) \
	FREE_ARRAY_PTR(lambdas) \
	thetas = NULL;

	/** MPI_Scatterv message:
	 *   [for each subproblem]
	 *   1 theta
	 *   2 lambda
	 */
	int size_of_sendbuf;     /**< MPI_Scatterv: size of send buffer pointer */
	double * sendbuf = NULL; /**< MPI_Scatterv: send buffer */
	int *    scounts = NULL; /**< MPI_Scatterv: send buffer size for each process */
	int *    sdispls = NULL; /**< MPI_Scatterv: send buffer displacement for each process*/

	/** MPI_Gatherv message:
	 *   [for each subproblem]
	 *   1 subproblem index
	 *   2 primal objective
	 *   3 dual objective
	 *   4 coupling column part of the solution
	 */
	int size_of_recvbuf;     /**< MPI_Gatherv: size of receive buffer pointer */
	double * recvbuf = NULL; /**< MPI_Gatherv: receive buffer */
	int *    rcounts = NULL; /**< MPI_Gatherv: receive buffer size for each process */
	int *    rdispls = NULL; /**< MPI_Gatherv: receive buffer displacement for each process*/

	const double * thetas  = NULL; /**< of master problem */
	double **      lambdas = NULL; /**< of master problem */

	BGN_TRY_CATCH

	int signal     = STO_STAT_MW_CONTINUE;   /**< signal to communicate with workers */
	int iter_count = 0;                      /**< iteration count */
	double ctime_start = CoinCpuTime();      /**< cputime of start */
	double wtime_start = CoinGetTimeOfDay(); /**< walltime of start */
	DecModel *   model   = master_->getModelPtr();   /**< pointer to model */
	StoMessage * message = master_->getMessagePtr(); /**< pointer to messsage object */

	/** allocate memory for buffer sizes and displacements */
	scounts = new int [comm_size_];
	sdispls = new int [comm_size_];
	rcounts = new int [comm_size_];
	rdispls = new int [comm_size_];

	/** initialize send buffer size and displacement,
	 * and calculate size of send buffer pointer */
	size_of_sendbuf = 0;
	for (int i = 0; i < comm_size_; ++i)
	{
		scounts[i] = 0;
		for (int j = 0; j < nsubprobs_[i]; ++j)
			scounts[i] += 1 + model->getNumSubproblemCouplingRows(subprob_indices_[subprob_displs_[i]+j]);
		sdispls[i] = i == 0 ? 0 : sdispls[i-1] + scounts[i-1];
		size_of_sendbuf += scounts[i];
	}

	/** initialize receive buffer size and displacement,
	 * and calculate size of receive buffer pointer */
	size_of_recvbuf = 0;
	rcounts[0] = 0;
	rdispls[0] = 0;
	for (int i = 1; i < comm_size_; ++i)
	{
		rcounts[i] = 0;
		for (int j = 0; j < nsubprobs_[i]; ++j)
			rcounts[i] += 3 + model->getNumSubproblemCouplingCols(i);
		rdispls[i] = rdispls[i-1] + rcounts[i-1];
		size_of_recvbuf += rcounts[i];
	}

	/** allocate memory for message buffers */
	sendbuf = new double [size_of_sendbuf];
	recvbuf = new double [size_of_recvbuf];

	/** allocate memory for lambdas */
	lambdas = new double * [model->getNumSubproblems()];

	/**
	 * This is the main loop to iteratively solve master problem.
	 */
	while (1)
	{
		/** receive message */
		MPI_Gatherv(NULL, 0, MPI_DOUBLE, recvbuf, rcounts, rdispls, MPI_DOUBLE, 0, comm_);

//		printf("master receive buffer:\n");
//		for (int i = 0; i < comm_size_; ++i)
//		{
//			printf("  rank %d:\n", i);
//			for (int j = 0; j < rcounts[i]; ++j)
//			{
//				if (j > 0 && j % 5 == 0) printf("\n");
//				printf("  [%5d] %+e", j, recvbuf[rdispls[i]+j]);
//			}
//			printf("\n");
//		}

		/** apply receive message */
		master_->nsubprobs_ = 0;
		for (int i = 0, j = 0, pos = 0; i < comm_size_; ++i)
		{
			message->print(3, "message count for rank %d: %d\n", i, rcounts[i]);
			for (int s = 0; s < nsubprobs_[i]; ++s)
			{
				master_->subindex_[j] = static_cast<int>(recvbuf[pos++]);
				master_->subprimobj_[j] = recvbuf[pos++];
				master_->subdualobj_[j] = recvbuf[pos++];
				CoinCopyN(recvbuf + pos, model->getNumSubproblemCouplingCols(master_->subindex_[j]), master_->subsolution_[j]);
				pos += model->getNumSubproblemCouplingCols(master_->subindex_[j]);
				message->print(5, "-> master, subprob %d primobj %+e\n", master_->subindex_[j], master_->subprimobj_[j]);
				j++;
			}
			master_->nsubprobs_ += nsubprobs_[i];
		}
		master_->worker_ = comm_size_ - 1;

		/** calculate dual objective */
		double dualobj = 0.0;
		for (int s = 0; s < master_->nsubprobs_; ++s)
			dualobj += master_->subprimobj_[s];
		master_->bestdualobj_ = dualobj > master_->bestdualobj_ ? dualobj : master_->bestdualobj_;
		message->print(2, "-> dual objective %e\n", dualobj);

		/** calculate absolute/relative gap */
		double absgap = fabs(master_->getPrimalObjective() - master_->getBestDualObjective());
		double relgap = absgap / (1.e-10 + fabs(master_->getPrimalObjective()));

		/** STOP with small gap */
		if (relgap < master_->getParPtr()->getDblParam("DD/STOP_TOL"))
		{
			signal = STO_STAT_MW_STOP;
			message->print(1, "STOP with gap tolerance %+e (%.2f%%).\n", absgap, relgap*100);
		}
		else if (iter_count > iteration_limit_)
		{
			signal = STO_STAT_MW_STOP;
			message->print(1, "STOP with iteration limit.\n");
		}
		else
		{
			/** update problem */
			master_->updateProblem();

			double tic = CoinGetTimeOfDay();

			/** solve problem */
			master_->solve();

			message->print(1, "Iteration %3d: Best primal %+e, Best dual %+e, Gap %+e (%.2f%%), Time elapsed %6.2f sec.\n",
					iter_count,
					master_->getBestPrimalObjective(), master_->getBestDualObjective(),
					absgap, relgap*100,
					CoinGetTimeOfDay() - wtime_start);

			/** increment iteration count */
			iter_count++;

			/** returns continue or stop signal */
			signal = master_->getStatus();
		}

		/** broadcast signal */
		MPI_Bcast(&signal, 1, MPI_INT, 0, comm_);

		/** signal to stop? */
		if (signal == STO_STAT_MW_STOP)
			break;

		/** retrieve master solution by part */
		double * master_primsol = const_cast<double*>(master_->getPrimalSolution());
		thetas  = master_primsol;
		for (int i = 0, j = model->getNumSubproblems(); i < model->getNumSubproblems(); ++i)
		{
			/** shallow copy */
			lambdas[i] = master_primsol + j;
			j += model->getNumSubproblemCouplingRows(i);
		}
		master_primsol = NULL;

		/** create send buffer */
		for (int i = 0, pos = 0; i < comm_size_; ++i)
		{
			for (int j = 0; j < nsubprobs_[i]; ++j)
			{
				int subprob_index = subprob_indices_[subprob_displs_[i]+j];
				sendbuf[pos++] = thetas[subprob_index];
				CoinCopyN(lambdas[subprob_index], model->getNumSubproblemCouplingRows(subprob_index), sendbuf + pos);
				pos += model->getNumSubproblemCouplingRows(subprob_index);
			}
		}

//		printf("master send buffer:\n");
//		for (int i = 0; i < comm_size_; ++i)
//		{
//			printf("  rank %d:\n", i);
//			for (int j = 0; j < scounts[i]; ++j)
//			{
//				if (j > 0 && j % 5 == 0) printf("\n");
//				printf("  [%5d] %e", j, sendbuf[sdispls[i]+j]);
//			}
//			printf("\n");
//		}

		/** send message */
		MPI_Scatterv(sendbuf, scounts, sdispls, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, 0, comm_);
	}

	/** set best dual objective */
	master_->bestdualobj_ = master_->primobj_;

//	printf("primsol_:\n");
//	for (int j = 0, k = 0; j < master_->getSiPtr()->getNumCols(); ++j)
//	{
//		if (fabs(master_->getPrimalSolution()[j]) < 1.0e-10) continue;
//		if (k > 0 && k % 5 == 0) printf("\n");
//		printf("  [%6d] %+e", j, master_->getPrimalSolution()[j]);
//		k++;
//	}
//	printf("\n");

	/** release shallow-copy of pointers */
	for (int i = 0; i < model->getNumSubproblems(); ++i)
		lambdas[i] = NULL;

	END_TRY_CATCH_RTN(FREE_MEMORY,STO_RTN_ERR)

	FREE_MEMORY

	return STO_RTN_OK;
#undef FREE_MEMORY
}

/** run master process */
STO_RTN_CODE DdMW::runMasterAsync()
{
#define FREE_MEMORY \
	FREE_ARRAY_PTR(sendbuf) \
	FREE_ARRAY_PTR(recvbuf) \
	FREE_ARRAY_PTR(lambdas)

	/** MPI_Send message:
	 *   [for each subproblem]
	 *   1 theta
	 *   2 lambda
	 */
	double * sendbuf = NULL; /**< MPI_Send: send buffer */
	int      scount  = 0;    /**< MPI_Send: send buffer size */

	/** MPI_Recv message:
	 *   [for each subproblem]
	 *   1 subproblem index
	 *   2 primal objective
	 *   3 dual objective
	 *   4 coupling column part of the solution
	 */
	double * recvbuf = NULL; /**< MPI_Recv: receive buffer */
	int      rcount  = 0;    /**< MPI_Recv: receive buffer size */

	const double * thetas  = NULL; /**< of master problem */
	double **      lambdas = NULL; /**< of master problem */

	int local_count; /**< local counter */

	BGN_TRY_CATCH

	MPI_Status status;
	int nworkers   = comm_size_ - 1;         /**< number of workers */
	int signal     = STO_STAT_MW_CONTINUE;   /**< signal to communicate with workers */
	int iter_count = 0;                      /**< iteration count */
	double ctime_start = CoinCpuTime();      /**< cputime of start */
	double wtime_start = CoinGetTimeOfDay(); /**< walltime of start */
	DecModel *   model   = master_->getModelPtr();   /**< pointer to model */
	StoMessage * message = master_->getMessagePtr(); /**< pointer to messsage object */

	/** maximum number of subproblems */
	int max_nsubprobs = 0;
	for (int i = 0; i < comm_size_; ++i)
		max_nsubprobs = max_nsubprobs < nsubprobs_[i] ? nsubprobs_[i] : max_nsubprobs;

	/** calculate send/receive buffer sizes */
	for (int i = 0; i < comm_size_; ++i)
	{
		int local_scount = 0, local_rcount = 0;
		for (int j = 0; j < nsubprobs_[i]; ++j)
		{
			local_scount += 1 + model->getNumSubproblemCouplingRows(subprob_indices_[subprob_displs_[i]+j]);
			local_rcount += 3 + model->getNumSubproblemCouplingCols(subprob_indices_[subprob_displs_[i]+j]);
		}
		scount = local_scount > scount ? local_scount : scount;
		rcount = local_rcount > rcount ? local_rcount : rcount;
	}

	/** allocate memory */
	sendbuf = new double [scount];
	recvbuf = new double [rcount];
	lambdas = new double * [model->getNumSubproblems()];

	while (signal == STO_STAT_MW_CONTINUE)
	{
		/** receive message */
		MPI_Recv(recvbuf, rcount, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, comm_, &status);

//		MPI_Get_count(&status, MPI_DOUBLE, &local_count);
//		printf("master receive buffer (%d):\n", local_count);
//		printf("  rank %d:\n", status.MPI_SOURCE);
//		for (int j = 0; j < local_count; ++j)
//		{
//			if (j > 0 && j % 5 == 0) printf("\n");
//			printf("  [%5d] %+e", j, recvbuf[j]);
//		}
//		printf("\n");

		/** signal to stop? */
		signal = status.MPI_TAG;
		if (signal == STO_STAT_MW_STOP)
			break;

		/** apply receive message */
		master_->worker_ = status.MPI_SOURCE;
		master_->nsubprobs_ = nsubprobs_[master_->worker_];
		for (int s = 0, pos = 0; s < master_->nsubprobs_; ++s)
		{
			master_->subindex_[s] = static_cast<int>(recvbuf[pos++]);
			master_->subprimobj_[s] = recvbuf[pos++];
			master_->subdualobj_[s] = recvbuf[pos++];
			CoinCopyN(recvbuf + pos, model->getNumSubproblemCouplingCols(master_->subindex_[s]), master_->subsolution_[s]);
			pos += model->getNumSubproblemCouplingCols(master_->subindex_[s]);
			message->print(5, "-> master received from rank %d: subprob %d primobj %+e\n",
					master_->worker_, master_->subindex_[s], master_->subprimobj_[s]);
		}

		/** update problem */
		master_->updateProblem();
		signal = master_->getStatus();

		if (signal == STO_STAT_MW_CONTINUE)
		{
			/** solve problem */
			double tic = CoinGetTimeOfDay();
			master_->solve();

			message->print(1, "Iteration %3d: Time elapsed %6.2f sec.\n",
					iter_count, CoinGetTimeOfDay() - wtime_start);

			/** increment iteration count */
			iter_count++;

			/** returns continue or stop signal */
			if (iter_count > iteration_limit_)
				signal = STO_STAT_MW_STOP;
			else
				signal = master_->getStatus();
		}

		/** send signal to stop */
		if (signal == STO_STAT_MW_STOP)
		{
			message->print(5, "-> send STOP signal to rank %d\n", master_->worker_);
			MPI_Send(NULL, 0, MPI_DOUBLE, master_->worker_, signal, comm_);
			break;
		}

		/** retrieve master solution by part */
		double * master_primsol = const_cast<double*>(master_->getPrimalSolution());
		thetas  = master_primsol;
		for (int i = 0, j = model->getNumSubproblems(); i < model->getNumSubproblems(); ++i)
		{
			/** shallow copy */
			lambdas[i] = master_primsol + j;
			j += model->getNumSubproblemCouplingRows(i);
		}
		master_primsol = NULL;

		/** create send buffer */
		int local_scount = 0;
		for (int s = 0; s < master_->nsubprobs_; ++s)
		{
			sendbuf[local_scount++] = thetas[master_->subindex_[s]];
			CoinCopyN(lambdas[master_->subindex_[s]], model->getNumSubproblemCouplingRows(master_->subindex_[s]), sendbuf + local_scount);
			local_scount += model->getNumSubproblemCouplingRows(master_->subindex_[s]);
		}

//		printf("master send buffer:\n");
//		printf("  rank %d:\n", master_->worker_);
//		for (int j = 0; j < local_scount; ++j)
//		{
//			if (j > 0 && j % 5 == 0) printf("\n");
//			printf("  [%6d] %+e", j, sendbuf[j]);
//		}
//		printf("\n");

		/** send message */
		MPI_Send(sendbuf, local_scount, MPI_DOUBLE, master_->worker_, signal, comm_);
	}

	/** sending stop signals one-by-one for asynchronous */
	nworkers--;

	/** sending stop signals to all the worker processes */
	double garbage = -1;
	while (nworkers > 0)
	{
		/** receive message */
		MPI_Recv(recvbuf, rcount, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, comm_, &status);
		/** send garbage with stop signal as a tag */
		MPI_Send(&garbage, 1, MPI_DOUBLE, status.MPI_SOURCE, signal, comm_);
		/** decrement number of workers */
		nworkers--;
	}

	/** set best dual objective */
	master_->bestdualobj_ = master_->primobj_;

	/** release shallow-copy of pointers */
	for (int i = 0; i < model->getNumSubproblems(); ++i)
		lambdas[i] = NULL;

	END_TRY_CATCH_RTN(FREE_MEMORY,STO_RTN_ERR)

	FREE_MEMORY

	return STO_RTN_OK;
#undef FREE_MEMORY
}

/** ruw worker processes */
STO_RTN_CODE DdMW::runWorker()
{
	if (comm_rank_ == 0)
		return STO_RTN_OK;

	DspParams * par = worker_->getParPtr();
	int narrprocidx = par->getIntPtrParamSize("ARR_PROC_IDX");
	double gaptol   = par->getDblParam("SCIP/GAP_TOL");

	for (int i = 0; i < narrprocidx; ++i)
		worker_->subprobs_[i]->setGapTol(100.0);

	if (sync_)
		return runWorkerSync();
	else
	{
		/** Initial round with rough tolerance */
		par->setDblParam("SCIP/GAP_TOL", 100.0);
		runWorkerSync();

		/** set the original tolerance */
		par->setDblParam("SCIP/GAP_TOL", gaptol);

		return runWorkerAsync();
	}
}

/** ruw worker processes */
STO_RTN_CODE DdMW::runWorkerSync()
{
#define FREE_MEMORY         \
	FREE_ARRAY_PTR(sendbuf) \
	FREE_ARRAY_PTR(recvbuf)

	/** MPI_Gatherv message:
	 *   [for each subproblem]
	 *   1 subproblem index
	 *   2 primal objective
	 *   3 dual objective
	 *   4 coupling column part of the solution
	 */
	double * sendbuf = NULL; /**< MPI_Gatherv: send buffer */
	int      scount  = 0;    /**< MPI_Gatherv: send buffer size */

	/** MPI_Scatterv message:
	 *   [for each subproblem]
	 *   1 theta
	 *   2 lambda
	 */
	double * recvbuf = NULL; /**< MPI_Scatterv: receive buffer */
	int      rcount  = 0;    /**< MPI_Scatterv: receive buffer size */

	BGN_TRY_CATCH

	int signal = STO_STAT_MW_CONTINUE;
	int narrprocidx  = worker_->getParPtr()->getIntPtrParamSize("ARR_PROC_IDX"); /**< number of subproblems */
	int * arrprocidx = worker_->getParPtr()->getIntPtrParam("ARR_PROC_IDX");     /**< subproblem indices */
	DecModel * model = worker_->getModelPtr();
	StoMessage * message = worker_->getMessagePtr(); /**< pointer to messsage object */

	/** calculate size of send buffer */
	for (int i = 0; i < narrprocidx; ++i)
		scount += 3 + model->getNumSubproblemCouplingCols(arrprocidx[i]);

	/** calculate size of receive buffer */
	for (int i = 0; i < narrprocidx; ++i)
		rcount += 1 + model->getNumSubproblemCouplingRows(arrprocidx[i]);

	/** allocate memory for message buffers */
	sendbuf = new double [scount];
	recvbuf = new double [rcount];

	/** loop until when the master signals stop */
	while (1)
	{
		/** Solve subproblems assigned to each process  */
		worker_->solve();

		/** TODO
		 * Add steps for generating Benders-type cuts and
		 * evaluating upper bound. If cuts are generated and
		 * exclude the current point, then we do not need to
		 * communicate with master.
		 *
		 * BendersWorker may be asynchronous.
		 */

		/** TODO send solution to Benders workers */

		/** TODO receive cuts and recourse function value */

		/** create send buffer */
		for (int s = 0, pos = 0; s < narrprocidx; ++s)
		{
			sendbuf[pos++] = static_cast<double>(worker_->subprobs_[s]->sind_);
			sendbuf[pos++] = worker_->subprobs_[s]->getPrimalBound();
			sendbuf[pos++] = worker_->subprobs_[s]->getDualBound();
			CoinCopyN(worker_->subprobs_[s]->si_->getSolution(), worker_->subprobs_[s]->ncols_coupling_, sendbuf + pos);
			pos += model->getNumSubproblemCouplingCols(worker_->subprobs_[s]->sind_);
			message->print(5, "-> worker %d, subprob %d primobj %+e dualobj %+e\n",
					comm_rank_, worker_->subprobs_[s]->sind_, worker_->subprobs_[s]->getPrimalBound(), worker_->subprobs_[s]->getDualBound());
		}

//		printf("Worker send message (%d):\n", scount);
//		for (int i = 0, j = 0; i < scount; ++i)
//		{
//			if (fabs(sendbuf[i]) < 1.0e-10) continue;
//			if (j > 0 && j % 5 == 0) printf("\n");
//			printf("  [%6d] %+e", i, sendbuf[i]);
//			j++;
//		}
//		printf("\n");

		/** send message to the master */
		MPI_Gatherv(sendbuf, scount, MPI_DOUBLE, NULL, NULL, NULL, MPI_DOUBLE, 0, comm_);

		/** receive message from the master */
		MPI_Bcast(&signal, 1, MPI_INT, 0, comm_);

		/** signal to stop? */
		if (signal == STO_STAT_MW_STOP)
			break;

		/** receive message from the master */
		MPI_Scatterv(NULL, NULL, NULL, MPI_DOUBLE, recvbuf, rcount, MPI_DOUBLE, 0, comm_);

//		printf("Worker received message (%d):\n", rcount);
//		for (int i = 0; i < rcount; ++i)
//		{
//			if (i > 0 && i % 5 == 0) printf("\n");
//			printf("  [%5d] %e", i, recvbuf[i]);
//		}
//		printf("\n");

		/** parse message */
		for (int s = 0, pos = 0; s < narrprocidx; ++s)
		{
			worker_->subprobs_[s]->theta_ = recvbuf[pos++];
			worker_->subprobs_[s]->updateProblem(recvbuf + pos);
			pos += model->getNumSubproblemCouplingRows(worker_->subprobs_[s]->sind_);
		}
	}

	/** release pointers */
	arrprocidx = NULL;
	model = NULL;

	END_TRY_CATCH_RTN(FREE_MEMORY,STO_RTN_ERR)

	FREE_MEMORY

	return STO_RTN_OK;
#undef FREE_MEMORY
}

/** ruw worker processes */
STO_RTN_CODE DdMW::runWorkerAsync()
{
#define FREE_MEMORY         \
	FREE_ARRAY_PTR(sendbuf) \
	FREE_ARRAY_PTR(recvbuf)

	/** MPI_Send message:
	 *   [for each subproblem]
	 *   1 subproblem index
	 *   2 primal objective
	 *   3 dual objective
	 *   4 coupling column part of the solution
	 */
	double * sendbuf = NULL; /**< MPI_Send: send buffer */
	int      scount  = 0;    /**< MPI_Send: send buffer size */

	/** MPI_Recv message:
	 *   [for each subproblem]
	 *   1 theta
	 *   2 lambda
	 */
	double * recvbuf = NULL; /**< MPI_Recv: receive buffer */
	int      rcount  = 0;    /**< MPI_Recv: receive buffer size */

	int local_count; /**< local counter */

	BGN_TRY_CATCH

	MPI_Status status;
	int signal       = STO_STAT_MW_CONTINUE; /**< signal to stop or continue */
	int narrprocidx  = worker_->getParPtr()->getIntPtrParamSize("ARR_PROC_IDX"); /**< number of subproblems */
	int * arrprocidx = worker_->getParPtr()->getIntPtrParam("ARR_PROC_IDX");     /**< subproblem indices */
	DecModel *   model   = worker_->getModelPtr();   /**< pointer to model object */
	StoMessage * message = worker_->getMessagePtr(); /**< pointer to message object */

	/** calculate size of send buffer */
	for (int i = 0; i < narrprocidx; ++i)
		scount += 3 + model->getNumSubproblemCouplingCols(arrprocidx[i]);

	/** calculate size of receive buffer */
	for (int i = 0; i < narrprocidx; ++i)
		rcount += 1 + model->getNumSubproblemCouplingRows(arrprocidx[i]);

	/** allocate memory for message buffers */
	sendbuf = new double [scount];
	recvbuf = new double [rcount];

	/** loop until when the master signals stop */
	while(1)
	{
		for (int i = 0; i < narrprocidx; ++i)
			worker_->subprobs_[i]->setGapTol(100.0);

		/** Solve subproblems assigned to each process  */
		worker_->solve();

		/** worker status */
		signal = worker_->getStatus();
		if (signal == STO_STAT_MW_STOP)
		{
			/** send stop signal to master */
			MPI_Send(NULL, 0, MPI_DOUBLE, 0, signal, comm_);
			break;
		}

		/** create send buffer */
		for (int s = 0, pos = 0; s < narrprocidx; ++s)
		{
			sendbuf[pos++] = static_cast<double>(worker_->subprobs_[s]->sind_);
			sendbuf[pos++] = worker_->subprobs_[s]->getPrimalBound();
			sendbuf[pos++] = worker_->subprobs_[s]->getDualBound();
			CoinCopyN(worker_->subprobs_[s]->si_->getSolution(), worker_->subprobs_[s]->ncols_coupling_, sendbuf + pos);
			pos += worker_->subprobs_[s]->ncols_coupling_;
			message->print(5, "-> worker %d, subprob %d primobj %+e dualobj %+e\n",
					comm_rank_, worker_->subprobs_[s]->sind_, worker_->subprobs_[s]->getPrimalBound(), worker_->subprobs_[s]->getDualBound());
		}

//		printf("Worker send message (%d):\n", scount);
//		for (int i = 0, j = 0; i < scount; ++i)
//		{
//			if (fabs(sendbuf[i]) < 1.0e-10) continue;
//			if (j > 0 && j % 5 == 0) printf("\n");
//			printf("  [%6d] %+e", i, sendbuf[i]);
//			j++;
//		}
//		printf("\n");

		/** send message to the master */
		MPI_Send(sendbuf, scount, MPI_DOUBLE, 0, signal, comm_);

		/** receive message from the master */
		MPI_Recv(recvbuf, rcount, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, comm_, &status);
		signal = status.MPI_TAG;

		if (signal == STO_STAT_MW_STOP)
			break;
		else
		{
//			MPI_Get_count(&status, MPI_DOUBLE, &local_count);
//			printf("Worker received message (%d):\n", local_count);
//			for (int i = 0, j = 0; i < local_count; ++i)
//			{
//				if (fabs(recvbuf[i]) < 1.0e-10) continue;
//				if (j > 0 && j % 5 == 0) printf("\n");
//				printf("  [%6d] %+e", i, recvbuf[i]);
//				j++;
//			}
//			printf("\n");
			/** parse message */
			for (int s = 0, pos = 0; s < narrprocidx; ++s)
			{
				worker_->subprobs_[s]->theta_ = recvbuf[pos++];
				worker_->subprobs_[s]->updateProblem(recvbuf + pos);
				pos += model->getNumSubproblemCouplingRows(worker_->subprobs_[s]->sind_);
			}
		}
	}

	/** release pointers */
	arrprocidx = NULL;
	model = NULL;

	END_TRY_CATCH_RTN(FREE_MEMORY,STO_RTN_ERR)

	FREE_MEMORY

	return STO_RTN_OK;
#undef FREE_MEMORY
}

STO_RTN_CODE DdMW::runBendersWorker()
{
	BGN_TRY_CATCH

	if (sync_)
	{
		/** Receive master solution */
//		MPI_Recv(solution, ncols, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, comm_, &status);
//		MPI_Bcast(solution, ncols, MPI_DOUBLE, 0, bdcomm_);

		/** TODO send cuts, if any, and recourse function value */
	}

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
}

/** finalize */
STO_RTN_CODE DdMW::finalize()
{
	BGN_TRY_CATCH

	/** collect results */
	/** collect statistics */

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
}

