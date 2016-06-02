/*
 * DdMWAsync.cpp
 *
 *  Created on: Mar 23, 2016
 *      Author: kibaekkim
 */

//#define DSP_DEBUG

#include "Solver/DualDecomp/DdMasterAtr.h"
#include "Solver/DualDecomp/DdMWAsync.h"

DdMWAsync::DdMWAsync(
		MPI_Comm          comm,   /**< MPI communicator */
		DdMaster *        master, /**< master problem */
		vector<DdWorker*> worker  /**< worker for finding lower bounds */):
DdMWPara(comm, master, worker) {}

DdMWAsync::~DdMWAsync() {}

DSP_RTN_CODE DdMWAsync::init()
{
	BGN_TRY_CATCH

	DdMWPara::init();
	sync_ = false;

	if (comm_rank_ == 0 && comm_size_ <= subcomm_size_)
		message_->print(1, "Upper bounding and the other extended features are disabled. (# of processors <= # of subprobs)\n");

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMWAsync::runMaster()
{
	if (comm_rank_ != 0)
		return DSP_RTN_OK;

	BGN_TRY_CATCH

	runMasterInit();
	runMasterCore();

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMWAsync::runWorker()
{
	if (comm_rank_ == 0)
		return DSP_RTN_OK;

	BGN_TRY_CATCH

	if (comm_color_ == 0)
	{
		DSPdebugMessage("Rank %d runs runWorkerInit() and runWorkerCore().\n", comm_rank_);
		runWorkerInit();
		runWorkerCore();
	}
	else
	{
		DSPdebugMessage("Rank %d runs runWorkerExt().\n", comm_rank_);
		runWorkerExt();
	}

	DSPdebugMessage("Rank %d finishied runWorker().\n", comm_rank_);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMWAsync::runMasterInit()
{
#define FREE_MEMORY \
	FREE_ARRAY_PTR(sendbuf) \
	FREE_ARRAY_PTR(scounts) \
	FREE_ARRAY_PTR(sdispls) \
	FREE_ARRAY_PTR(recvbuf) \
	FREE_ARRAY_PTR(rcounts) \
	FREE_ARRAY_PTR(rdispls) \
	FREE_ARRAY_PTR(lambdas) \
	FREE_2D_ARRAY_PTR(subcomm_size_-1,subindex)   \
	FREE_2D_ARRAY_PTR(subcomm_size_-1,subprimobj) \
	FREE_2D_ARRAY_PTR(subcomm_size_-1,subdualobj) \
	if (subsolution) {                                         \
		for (int i = 0; i < subcomm_size_ - 1; ++i) {          \
			FREE_2D_ARRAY_PTR(nsubprobs_[i+1], subsolution[i]) \
		}                                                      \
		delete [] subsolution;                                 \
		subsolution = NULL;                                    \
	}                                                          \
	FREE_ARRAY_PTR(number_of_elements) \
	FREE_ARRAY_PTR(indices)            \
	FREE_ARRAY_PTR(elements)

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

	int ** subindex = NULL;
	double ** subprimobj = NULL;
	double ** subdualobj = NULL;
	double *** subsolution = NULL;

	/** to send solutions to workers */
	MPI_Status status;
	MPI_Request request;
	int *    number_of_elements = NULL; /**< number of nonzero elements in solution vector */
	int *    indices            = NULL; /**< indices of solution vectors */
	double * elements           = NULL; /**< elements of solution vectors */

	BGN_TRY_CATCH

	DdMasterAtr * master = dynamic_cast<DdMasterAtr*>(master_);

	/** allocate memory for buffer sizes and displacements */
	scounts = new int [subcomm_size_];
	sdispls = new int [subcomm_size_];
	rcounts = new int [subcomm_size_];
	rdispls = new int [subcomm_size_];

	/** initialize send buffer size and displacement,
	 * and calculate size of send buffer pointer */
	size_of_sendbuf = 0;
	for (int i = 0; i < subcomm_size_; ++i)
	{
		scounts[i] = 0;
		for (int j = 0; j < nsubprobs_[i]; ++j)
			scounts[i] += 1 + model_->getNumSubproblemCouplingRows(subprob_indices_[subprob_displs_[i]+j]);
		sdispls[i] = i == 0 ? 0 : sdispls[i-1] + scounts[i-1];
		size_of_sendbuf += scounts[i];
	}

	/** initialize receive buffer size and displacement,
	 * and calculate size of receive buffer pointer */
	size_of_recvbuf = 0;
	rcounts[0] = 0;
	rdispls[0] = 0;
	for (int i = 1; i < subcomm_size_; ++i)
	{
		rcounts[i] = 0;
		for (int j = 0; j < nsubprobs_[i]; ++j)
			rcounts[i] += 3 + model_->getNumSubproblemCouplingCols(i);
		rdispls[i] = rdispls[i-1] + rcounts[i-1];
		size_of_recvbuf += rcounts[i];
	}

	/** allocate memory for message buffers */
	sendbuf = new double [size_of_sendbuf];
	recvbuf = new double [size_of_recvbuf];

	/** allocate memory for lambdas */
	lambdas = new double * [model_->getNumSubproblems()];

	/** allocate memory for subproblem results */
	subindex    = new int * [subcomm_size_ - 1];
	subprimobj  = new double * [subcomm_size_ - 1];
	subdualobj  = new double * [subcomm_size_ - 1];
	subsolution = new double ** [subcomm_size_ - 1];
	for (int i = 0; i < subcomm_size_ - 1; ++i)
	{
		subindex[i]    = new int [nsubprobs_[i+1]];
		subprimobj[i]  = new double [nsubprobs_[i+1]];
		subdualobj[i]  = new double [nsubprobs_[i+1]];
		subsolution[i] = new double * [nsubprobs_[i+1]];
		for (int s1 = 0, s2 = subprob_displs_[i+1]; s1 < nsubprobs_[i+1]; ++s1)
		{
			subsolution[i][s1] = new double [model_->getNumSubproblemCouplingCols(subprob_indices_[s2])];
			s2++;
		}
	}

	/** to send solutions to workers */
	int maxnum_of_solutions = subcomm_size_; /**< maximum number of solutions */
	int maxnum_of_indices   = subcomm_size_ * model_->getNumSubproblemCouplingCols(0); /**< maximum length of indices (elements) */
	number_of_elements = new int [maxnum_of_solutions];
	indices            = new int [maxnum_of_indices];
	elements           = new double [maxnum_of_indices];

	/** receive message */
	MPI_Gatherv(NULL, 0, MPI_DOUBLE, recvbuf, rcounts, rdispls, MPI_DOUBLE, 0, subcomm_);

//	printf("master receive buffer:\n");
//	for (int i = 0; i < subcomm_size_; ++i)
//	{
//		printf("  rank %d:\n", i);
//		StoMessage::printArray(rcounts[i], recvbuf + rdispls[i]);
//	}

	/** apply receive message */
	double dualobj = 0.0;
	for (int i = 1, pos = 0; i < subcomm_size_; ++i)
	{
		DSPdebugMessage("message count for rank %d: %d\n", i, rcounts[i]);
		for (int s = 0; s < nsubprobs_[i]; ++s)
		{
			subindex[i-1][s] = static_cast<int>(recvbuf[pos++]);
			subprimobj[i-1][s] = recvbuf[pos++];
			subdualobj[i-1][s] = recvbuf[pos++];
			CoinCopyN(recvbuf + pos,
					model_->getNumSubproblemCouplingCols(subindex[i-1][s]),
					subsolution[i-1][s]);
			pos += model_->getNumSubproblemCouplingCols(subindex[i-1][s]);
			DSPdebugMessage("-> master, subprob %d primobj %+e\n",
					subindex[i-1][s], subprimobj[i-1][s]);
			dualobj += subprimobj[i-1][s];
		}
		master->worker_.push_back(i);
		master->nsubprobs_.push_back(nsubprobs_[i]);
		master->subindex_.push_back(subindex[i-1]);
		master->subprimobj_.push_back(subprimobj[i-1]);
		master->subdualobj_.push_back(subdualobj[i-1]);
		master->subsolution_.push_back(subsolution[i-1]);
	}

	int num_ext_procs = comm_size_ - subcomm_size_;
	if (num_ext_procs < 0) num_ext_procs = 0;

	if (num_ext_procs > 0)
	{
		/** number of solutions in the pool */
		int nSolutionsPool = ubSolutions_.size();

		/** send coupling solution for upper bounds and cuts */
		sendCouplingSolutions(maxnum_of_solutions, maxnum_of_indices, number_of_elements, indices, elements, &request);

		int nSolutionsAdded = ubSolutions_.size() - nSolutionsPool;
		message_->print(3, "-> %d solutions added to the pool (%d)\n", nSolutionsAdded, nSolutionsPool + nSolutionsAdded);
	}

	/** calculate dual objective */
	master_->bestdualobj_ = dualobj > master_->bestdualobj_ ? dualobj : master_->bestdualobj_;

	/** update problem */
	master_->updateProblem();

	/** solve problem */
	master_->solve();

	/** retrieve master solution by part */
	double * master_primsol = const_cast<double*>(master_->getPrimalSolution());
	thetas  = master_primsol;
	for (int i = 0, j = model_->getNumSubproblems(); i < model_->getNumSubproblems(); ++i)
	{
		/** shallow copy */
		lambdas[i] = master_primsol + j;
		j += model_->getNumSubproblemCouplingRows(i);
	}
	master_primsol = NULL;

	/** create send buffer */
	for (int i = 0, pos = 0; i < subcomm_size_; ++i)
	{
		for (int j = 0; j < nsubprobs_[i]; ++j)
		{
			int subprob_index = subprob_indices_[subprob_displs_[i]+j];
			sendbuf[pos++] = thetas[subprob_index];
			CoinCopyN(lambdas[subprob_index], model_->getNumSubproblemCouplingRows(subprob_index), sendbuf + pos);
			pos += model_->getNumSubproblemCouplingRows(subprob_index);
		}
	}

//	printf("master send buffer:\n");
//	for (int i = 0; i < subcomm_size_; ++i)
//	{
//		printf("  rank %d:\n", i);
//		StoMessage::printArray(scounts[i], sendbuf + sdispls[i]);
//	}

	/** send message */
	MPI_Scatterv(sendbuf, scounts, sdispls, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, 0, subcomm_);

	/** release shallow-copy of pointers */
	for (int i = 0; i < model_->getNumSubproblems(); ++i)
		lambdas[i] = NULL;

	/** make sure if the solutions were safely buffered. */
	MPI_Wait(&request, &status);

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

DSP_RTN_CODE DdMWAsync::runMasterCore()
{
#define FREE_MEMORY \
	FREE_ARRAY_PTR(sendbuf) \
	FREE_ARRAY_PTR(recvbuf) \
	FREE_ARRAY_PTR(lambdas) \
	FREE_2D_ARRAY_PTR(subcomm_size_-1,subindex)   \
	FREE_2D_ARRAY_PTR(subcomm_size_-1,subprimobj) \
	FREE_2D_ARRAY_PTR(subcomm_size_-1,subdualobj) \
	if (subsolution) {                                                  \
		for (int i = 0; i < subcomm_size_ - 1; ++i) {                   \
			FREE_2D_ARRAY_PTR(master_->maxnumsubprobs_, subsolution[i]) \
		}                                                               \
		delete [] subsolution;                                          \
		subsolution = NULL;                                             \
	}                                                                   \
	FREE_ARRAY_PTR(primobjs)           \
	FREE_ARRAY_PTR(number_of_elements) \
	FREE_ARRAY_PTR(indices)            \
	FREE_ARRAY_PTR(elements)           \
	FREE_ARRAY_PTR(subprimobjFilled)

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

	int ** subindex = NULL;
	double ** subprimobj = NULL;
	double ** subdualobj = NULL;
	double *** subsolution = NULL;

	int nSolutionsPool;  /**< number of solutions in the solution pool */
	int nSolutionsAdded; /**< number of solutions added to the pool */
	double * primobjs     = NULL; /**< primal objectives from subproblems */

	char itercode;
	int recv_message;

	int local_count; /**< local counter */

	/** to send solutions to workers */
	int *    number_of_elements = NULL; /**< number of nonzero elements in solution vector */
	int *    indices            = NULL; /**< indices of solution vectors */
	double * elements           = NULL; /**< elements of solution vectors */

	bool * subprimobjFilled = NULL; /**< indicate if subproblem objective value is filled for each worker */

	BGN_TRY_CATCH

	MPI_Status status;
	MPI_Request request = MPI_REQUEST_NULL;
	int nworkers   = subcomm_size_ - 1;      /**< number of workers */
	int signal     = DSP_STAT_MW_CONTINUE;   /**< signal to communicate with workers */
	int iter_count = 0;                      /**< iteration count */
	double ctime_start = CoinCpuTime();      /**< cputime of start */
	double wtime_start = CoinGetTimeOfDay(); /**< walltime of start */
	DdMasterAtr * master = dynamic_cast<DdMasterAtr*>(master_);

	/** maximum number of subproblems */
	int max_nsubprobs = 0;
	for (int i = 0; i < subcomm_size_; ++i)
		max_nsubprobs = max_nsubprobs < nsubprobs_[i] ? nsubprobs_[i] : max_nsubprobs;

	/** calculate send/receive buffer sizes */
	for (int i = 0; i < subcomm_size_; ++i)
	{
		int local_scount = 0, local_rcount = 0;
		for (int j = 0; j < nsubprobs_[i]; ++j)
		{
			local_scount += 1 + model_->getNumSubproblemCouplingRows(subprob_indices_[subprob_displs_[i]+j]);
			local_rcount += 3 + model_->getNumSubproblemCouplingCols(subprob_indices_[subprob_displs_[i]+j]);
		}
		scount = local_scount > scount ? local_scount : scount;
		rcount = local_rcount > rcount ? local_rcount : rcount;
	}

	/** allocate memory */
	sendbuf = new double [scount];
	recvbuf = new double [rcount];
	lambdas = new double * [model_->getNumSubproblems()];

	/** allocate memory for subproblem results */
	subindex    = new int * [subcomm_size_ - 1];
	subprimobj  = new double * [subcomm_size_ - 1];
	subdualobj  = new double * [subcomm_size_ - 1];
	subsolution = new double ** [subcomm_size_ - 1];
	for (int i = 0; i < subcomm_size_ - 1; ++i)
	{
		subindex[i]    = new int [master_->maxnumsubprobs_];
		subprimobj[i]  = new double [master_->maxnumsubprobs_];
		subdualobj[i]  = new double [master_->maxnumsubprobs_];
		subsolution[i] = new double * [master_->maxnumsubprobs_];
		for (int s = 0; s < master_->maxnumsubprobs_; ++s)
			subsolution[i][s] = new double [master_->getModelPtr()->getNumSubproblemCouplingCols(s)];
	}

	/** allocate memory */
	int primobjs_capacity = model_->getNumSubproblems();
	primobjs     = new double [primobjs_capacity];

	int num_ext_procs = comm_size_ - subcomm_size_;
	if (num_ext_procs < 0) num_ext_procs = 0;
	DSPdebugMessage("Number of processors for DdWorkerExt(): %d\n", num_ext_procs);

	/** to send solutions to workers */
	int maxnum_of_solutions = subcomm_size_; /**< maximum number of solutions */
	int maxnum_of_indices   = subcomm_size_ * model_->getNumSubproblemCouplingCols(0); /**< maximum length of indices (elements) */
	number_of_elements = new int [maxnum_of_solutions];
	indices            = new int [maxnum_of_indices];
	elements           = new double [maxnum_of_indices];

	bool exactEval = false;
	bool checkpointEnabled = false; /**< indicate if all the subproblems are given the same dual values. */
	int nsubprimobjFilled = subcomm_size_ - 1;
	subprimobjFilled = new bool [subcomm_size_ - 1];
	CoinFillN(subprimobjFilled, subcomm_size_ - 1, true);
	double dualobj = 0.0;

	/**
	 * PRINT DISPLAY
	 *
	 * 0123456789012345678901234567890123456789
	 * iter      curobj     primobj dualobj gap time
	 *  %4d  %+10e  %+10e  %+10e  %6.2f  %6.1f
	 */
	message_->print(1, "  %4s  %13s  %13s  %13s  %6s  %6s\n",
			"iter", "curobj", "primobj", "dualobj", "gap(%)", "times");

	while (signal != DSP_STAT_MW_STOP)
	{
		itercode = ' ';
		recv_message = 1;

		/** clear queue */
		master->worker_.clear();
		master->nsubprobs_.clear();
		master->subindex_.clear();
		master->subprimobj_.clear();
		master->subdualobj_.clear();
		master->subsolution_.clear();

//		printf("################# time stamp 1: %.2f\n", CoinGetTimeOfDay() - wtime_start);
		while (recv_message)
		{
			/** receive message */
			MPI_Recv(recvbuf, rcount, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, subcomm_, &status);

//			DSPdebug(MPI_Get_count(&status, MPI_DOUBLE, &local_count));
//			DSPdebugMessage("master receive buffer (%d):\n", local_count);
//			DSPdebugMessage("  rank %d:\n", status.MPI_SOURCE);
//			DSPdebug(StoMessage::printArray(local_count, recvbuf));

			/** signal to stop? */
			signal = status.MPI_TAG;
			if (signal == DSP_STAT_MW_STOP)
				break;

			/** apply receive message */
			master->worker_.push_back(status.MPI_SOURCE);
			master->nsubprobs_.push_back(nsubprobs_[status.MPI_SOURCE]);
			for (int s = 0, pos = 0; s < nsubprobs_[status.MPI_SOURCE]; ++s)
			{
				subindex[status.MPI_SOURCE-1][s] = static_cast<int>(recvbuf[pos++]);
				subprimobj[status.MPI_SOURCE-1][s] = recvbuf[pos++];
				subdualobj[status.MPI_SOURCE-1][s] = recvbuf[pos++];
				CoinCopyN(recvbuf + pos,
						model_->getNumSubproblemCouplingCols(subindex[status.MPI_SOURCE-1][s]),
						subsolution[status.MPI_SOURCE-1][s]);
				pos += model_->getNumSubproblemCouplingCols(subindex[status.MPI_SOURCE-1][s]);
				DSPdebugMessage("master received from rank %d: subprob %d primobj %+e\n",
						status.MPI_SOURCE, subindex[status.MPI_SOURCE-1][s], subprimobj[status.MPI_SOURCE-1][s]);
			}
			master->subindex_.push_back(subindex[status.MPI_SOURCE-1]);
			master->subprimobj_.push_back(subprimobj[status.MPI_SOURCE-1]);
			master->subdualobj_.push_back(subdualobj[status.MPI_SOURCE-1]);
			master->subsolution_.push_back(subsolution[status.MPI_SOURCE-1]);

			DSPdebugMessage("worker size %lu, cutoff %.2f\n", master->worker_.size(), (subcomm_size_ - 1) / 2.0);
//			if (checkpointEnabled == true || master->worker_.size() < (subcomm_size_ - 1) / 2.0)
			{
				MPI_Status probe_status;
				MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, subcomm_, &recv_message, &probe_status);
				if (recv_message)
				{
					MPI_Get_count(&probe_status, MPI_DOUBLE, &local_count);
					DSPdebugMessage("Iprobe: source %d count %d recv_message %d\n", probe_status.MPI_SOURCE, local_count, recv_message);
				}
			}

			if (master->worker_.size() == subcomm_size_ - 1)
				break;
		}
//		printf("################# time stamp 2: %.2f\n", CoinGetTimeOfDay() - wtime_start);

		DSPdebugMessage("Number of worker messages: %lu\n", master->worker_.size());
		if (checkpointEnabled)
		{
			for (unsigned i = 0; i < master->worker_.size(); ++i)
			{
				if (subprimobjFilled[master->worker_[i]-1] == false)
				{
					for (int s = 0; s < nsubprobs_[master->worker_[i]]; ++s)
					{
						DSPdebugMessage("worker %d dualobj %+e\n", master->worker_[i], master->subdualobj_[master->worker_[i]-1][s]);
						dualobj += master->subdualobj_[master->worker_[i]-1][s];
					}
					subprimobjFilled[master->worker_[i]-1] = true;
					nsubprimobjFilled++;
				}
			}
			if (nsubprimobjFilled == subcomm_size_ - 1)
			{
				DSPdebugMessage("Found new dual objective %+e, disabled check point\n", dualobj);
				checkpointEnabled = false;
				if (dualobj > master_->bestdualobj_)
				{
					master_->bestdualobj_ = dualobj;
					itercode = '*';
				}
			}
		}
		else if (master->worker_.size() == subcomm_size_ - 1)
		{
			DSPdebugMessage("Enabled check point\n");
			checkpointEnabled = true;
			exactEval = true;
			CoinFillN(subprimobjFilled, subcomm_size_ - 1, false);
			nsubprimobjFilled = 0;
			dualobj = 0.0;
		}

		/**
		 * send a set of solutions; and receive upper bounds
		 */
		if (num_ext_procs > 0)
		{
			/** number of solutions in the pool */
			nSolutionsPool = ubSolutions_.size();

			/** make sure if the solutions were safely buffered. */
			if (request != MPI_REQUEST_NULL)
				MPI_Wait(&request, &status);

			/** send coupling solution for upper bounds and cuts */
			sendCouplingSolutions(maxnum_of_solutions, maxnum_of_indices,
					number_of_elements, indices, elements, &request);

			nSolutionsAdded = ubSolutions_.size() - nSolutionsPool;
			message_->print(3, "-> %d solutions added to the pool (%d)\n", nSolutionsAdded, nSolutionsPool + nSolutionsAdded);

			/** is there a message to receive? */
			MPI_Iprobe(comm_root_key_, MPI_ANY_TAG, comm_, &recv_message, &status);
			DSPdebugMessage("MPI_Iprobe returns recv_message %d.\n" ,recv_message);

			/** receive upper bounds and cuts */
			double bestprimobj = master_->bestprimobj_;
			while (recv_message)
			{
				int num_primobjs;
				MPI_Get_count(&status, MPI_DOUBLE, &num_primobjs);

				/** realloc memeory */
				if (num_primobjs > primobjs_capacity)
				{
					primobjs = (double*) realloc(primobjs, num_primobjs * sizeof(double));
					primobjs_capacity = num_primobjs;
				}

				MPI_Recv(primobjs, num_primobjs, MPI_DOUBLE, comm_root_key_, MPI_ANY_TAG, comm_, &status);

				/** calculate primal objective */
				for (int i = 0; i < num_primobjs; ++i)
				{
					DSPdebugMessage("solution %d: primal objective %+e\n", i, primobjs[i]);
					master_->bestprimobj_ = primobjs[i] < master_->bestprimobj_ ? primobjs[i] : master_->bestprimobj_;
				}

				/** is there a message to receive? */
				MPI_Iprobe(comm_root_key_, MPI_ANY_TAG, comm_, &recv_message, &status);
			}

			if (bestprimobj > master_->bestprimobj_)
				itercode = '*';
		}
//		printf("################# time stamp 3: %.2f\n", CoinGetTimeOfDay() - wtime_start);

		/** calculate absolute/relative gap */
		double absgap = fabs(master_->getBestPrimalObjective() - master_->getBestDualObjective());
		double relgap = absgap / (1.e-10 + fabs(master_->getBestPrimalObjective()));

		/** update problem */
		master->updateProblem();

		/** solve problem */
		double master_time = CoinGetTimeOfDay();
		master->solve();

		message_->print(3, "-> master solution time %.1f sec.\n", CoinGetTimeOfDay() - master_time);

		DSPdebugMessage("iteration %d: master status %d\n", iter_count, master->getStatus());

		message_->print(1, " %c%4d  %+10e", itercode, iter_count, master_->getPrimalObjective());
		if (master_->getBestPrimalObjective() < 1.0e+20)
			message_->print(1, "  %+10e", master_->getBestPrimalObjective());
		else
			message_->print(1, "  %+13s", "inf");
		message_->print(1, "  %+10e  %6.2f  %6.1f\n",
				master_->getBestDualObjective(), relgap*100, CoinGetTimeOfDay() - wtime_start);

		/** increment iteration count */
		iter_count++;

		/** returns continue or stop signal */
		if (iter_count > master_->getParPtr()->getIntParam("ITER_LIM"))
			signal = DSP_STAT_MW_STOP;
		else
			signal = master->terminationTest();

		/** send signal to stop */
		if (signal == DSP_STAT_MW_STOP)
		{
			for (unsigned i = 0; i < master->worker_.size(); ++i)
			{
				message_->print(5, "-> send STOP signal to rank %d (%d)\n", master->worker_[i], nworkers);
				MPI_Send(NULL, 0, MPI_DOUBLE, master->worker_[i], signal, subcomm_);
				nworkers--;
			}
			if (num_ext_procs > 0)
			{
				int czero = 0;
				message_->print(5, "-> send STOP signal to rank %d\n", comm_root_key_);
				MPI_Send(&czero, 1, MPI_INT, comm_root_key_, DSP_STAT_MW_STOP, comm_);
			}
			break;
		}

		/** retrieve master solution by part */
		double * master_primsol = const_cast<double*>(master_->getPrimalSolution());
		thetas  = master_primsol;
		for (int i = 0, j = model_->getNumSubproblems(); i < model_->getNumSubproblems(); ++i)
		{
			/** shallow copy */
			lambdas[i] = master_primsol + j;
			j += model_->getNumSubproblemCouplingRows(i);
		}
		master_primsol = NULL;

		/** create and send message buffer */

		if (exactEval == true && signal == DSP_STAT_MW_CONTINUE)
		{
			signal = DSP_STAT_MW_EXACT;
			exactEval = false;
		}
		for (unsigned i = 0; i < master->worker_.size(); ++i)
		{
			int local_scount = 0;
			for (unsigned s = 0; s < master->nsubprobs_[i]; ++s)
			{
				sendbuf[local_scount++] = thetas[master->subindex_[i][s]];
				CoinCopyN(lambdas[master->subindex_[i][s]], model_->getNumSubproblemCouplingRows(master->subindex_[i][s]), sendbuf + local_scount);
				local_scount += model_->getNumSubproblemCouplingRows(master->subindex_[i][s]);
			}

//			DSPdebugMessage("master send buffer:\n");
//			DSPdebug(StoMessage::printArray(local_scount, sendbuf));
			MPI_Send(sendbuf, local_scount, MPI_DOUBLE, master->worker_[i], signal, subcomm_);
		}
	}

	/** sending stop signals to all the worker processes */
	double garbage = -1;
	while (nworkers > 0)
	{
		/** receive message */
		MPI_Recv(recvbuf, rcount, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, subcomm_, &status);
		/** send garbage with stop signal as a tag */
		message_->print(5, "-> send STOP signal to rank %d (%d)\n", status.MPI_SOURCE, nworkers);
		MPI_Send(&garbage, 1, MPI_DOUBLE, status.MPI_SOURCE, DSP_STAT_MW_STOP, subcomm_);
		/** decrement number of workers */
		nworkers--;
	}
	DSPdebugMessage("finishing runDdMAsterCore().\n");

	/** clear queue */
	master->worker_.clear();
	master->nsubprobs_.clear();
	master->subindex_.clear();
	master->subprimobj_.clear();
	master->subdualobj_.clear();
	master->subsolution_.clear();

	/** set best dual objective */
	master_->bestdualobj_ = master_->primobj_;

	/** release shallow-copy of pointers */
	for (int i = 0; i < model_->getNumSubproblems(); ++i)
		lambdas[i] = NULL;

	/** make sure if the solutions were safely buffered. */
	if (request != MPI_REQUEST_NULL)
		MPI_Wait(&request, &status);

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

DSP_RTN_CODE DdMWAsync::runWorkerInit()
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

	assert(worker_[0]->getType()==DdWorker::LB);
	DdWorkerLB * workerlb = dynamic_cast<DdWorkerLB*>(worker_[0]);

	int narrprocidx  = workerlb->getParPtr()->getIntPtrParamSize("ARR_PROC_IDX"); /**< number of subproblems */
	int * arrprocidx = workerlb->getParPtr()->getIntPtrParam("ARR_PROC_IDX");     /**< subproblem indices */

	/** calculate size of send buffer */
	for (int i = 0; i < narrprocidx; ++i)
		scount += 3 + model_->getNumSubproblemCouplingCols(arrprocidx[i]);

	/** calculate size of receive buffer */
	for (int i = 0; i < narrprocidx; ++i)
		rcount += 1 + model_->getNumSubproblemCouplingRows(arrprocidx[i]);

	/** allocate memory for message buffers */
	sendbuf = new double [scount];
	recvbuf = new double [rcount];

	/** Solve subproblems assigned to each process  */
	workerlb->solve();

	/** create send buffer */
	for (int s = 0, pos = 0; s < narrprocidx; ++s)
	{
		sendbuf[pos++] = static_cast<double>(workerlb->subprobs_[s]->sind_);
		sendbuf[pos++] = workerlb->subprobs_[s]->getPrimalBound();
		sendbuf[pos++] = workerlb->subprobs_[s]->getDualBound();
		CoinCopyN(workerlb->subprobs_[s]->si_->getSolution(), workerlb->subprobs_[s]->ncols_coupling_, sendbuf + pos);
		pos += model_->getNumSubproblemCouplingCols(workerlb->subprobs_[s]->sind_);
		message_->print(5, "MW -> worker %d, subprob %d primobj %+e dualobj %+e\n",
				comm_rank_, workerlb->subprobs_[s]->sind_, workerlb->subprobs_[s]->getPrimalBound(), workerlb->subprobs_[s]->getDualBound());
	}

//	DSPdebugMessage("Worker send message (%d):\n", scount);
//	DSPdebug(StoMessage::printArray(scount, sendbuf));

	/** send message to the master */
	MPI_Gatherv(sendbuf, scount, MPI_DOUBLE, NULL, NULL, NULL, MPI_DOUBLE, 0, subcomm_);

	/** receive message from the master */
	MPI_Scatterv(NULL, NULL, NULL, MPI_DOUBLE, recvbuf, rcount, MPI_DOUBLE, 0, subcomm_);

//	DSPdebugMessage("Worker received message (%d):\n", rcount);
//	DSPdebug(StoMessage::printArray(rcount, recvbuf));

	/** parse message */
	for (int s = 0, pos = 0; s < narrprocidx; ++s)
	{
		workerlb->subprobs_[s]->theta_ = recvbuf[pos++];
		workerlb->subprobs_[s]->updateProblem(recvbuf + pos);
		pos += model_->getNumSubproblemCouplingRows(workerlb->subprobs_[s]->sind_);
	}

	/** release pointers */
	arrprocidx = NULL;

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

DSP_RTN_CODE DdMWAsync::runWorkerCore()
{
	if (worker_.size() == 0)
		return DSP_RTN_OK;

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

	assert(worker_[0]->getType()==DdWorker::LB);
	DdWorkerLB * workerlb = dynamic_cast<DdWorkerLB*>(worker_[0]);

	MPI_Status status;
	int signal       = DSP_STAT_MW_CONTINUE; /**< signal to stop or continue */
	int narrprocidx  = workerlb->getParPtr()->getIntPtrParamSize("ARR_PROC_IDX"); /**< number of subproblems */
	int * arrprocidx = workerlb->getParPtr()->getIntPtrParam("ARR_PROC_IDX");     /**< subproblem indices */

	/** calculate size of send buffer */
	for (int i = 0; i < narrprocidx; ++i)
		scount += 3 + model_->getNumSubproblemCouplingCols(arrprocidx[i]);

	/** calculate size of receive buffer */
	for (int i = 0; i < narrprocidx; ++i)
		rcount += 1 + model_->getNumSubproblemCouplingRows(arrprocidx[i]);

	/** allocate memory for message buffers */
	sendbuf = new double [scount];
	recvbuf = new double [rcount];

	/** loop until when the master signals stop */
	while(1)
	{
		double gaptol = 100.0;
		if (signal == DSP_STAT_MW_EXACT)
			gaptol = par_->getDblParam("DD/STOP_TOL");

		/** set gap tolerance */
		for (int i = 0; i < narrprocidx; ++i)
			workerlb->subprobs_[i]->setGapTol(gaptol);

		/** Solve subproblems assigned to each process  */
		workerlb->solve();

		/** worker status */
		signal = workerlb->getStatus();
		if (signal == DSP_STAT_MW_STOP)
		{
			/** send stop signal to master */
			MPI_Send(NULL, 0, MPI_DOUBLE, 0, signal, subcomm_);
			break;
		}

		/** create send buffer */
		for (int s = 0, pos = 0; s < narrprocidx; ++s)
		{
			sendbuf[pos++] = static_cast<double>(workerlb->subprobs_[s]->sind_);
			sendbuf[pos++] = workerlb->subprobs_[s]->getPrimalBound();
			sendbuf[pos++] = workerlb->subprobs_[s]->getDualBound();
			CoinCopyN(workerlb->subprobs_[s]->si_->getSolution(), workerlb->subprobs_[s]->ncols_coupling_, sendbuf + pos);
			pos += workerlb->subprobs_[s]->ncols_coupling_;
			DSPdebugMessage("worker %d, subprob %d primobj %+e dualobj %+e\n",
					comm_rank_, workerlb->subprobs_[s]->sind_, workerlb->subprobs_[s]->getPrimalBound(), workerlb->subprobs_[s]->getDualBound());
		}

//		DSPdebugMessage("Worker send message (%d):\n", scount);
//		DSPdebug(StoMessage::printArray(scount, sendbuf));

		/** send message to the master */
		MPI_Send(sendbuf, scount, MPI_DOUBLE, 0, signal, subcomm_);

		/** receive message from the master */
		MPI_Recv(recvbuf, rcount, MPI_DOUBLE, 0, MPI_ANY_TAG, subcomm_, &status);
		signal = status.MPI_TAG;

		if (signal == DSP_STAT_MW_STOP)
			break;
		else
		{
//			DSPdebug(MPI_Get_count(&status, MPI_DOUBLE, &local_count));
//			DSPdebugMessage("Worker received message (%d):\n", local_count);
//			DSPdebug(StoMessage::printArray(local_count, recvbuf));

			/** parse message */
			for (int s = 0, pos = 0; s < narrprocidx; ++s)
			{
				workerlb->subprobs_[s]->theta_ = recvbuf[pos++];
				workerlb->subprobs_[s]->updateProblem(recvbuf + pos);
				pos += model_->getNumSubproblemCouplingRows(workerlb->subprobs_[s]->sind_);
			}
		}
	}

	/** release pointers */
	arrprocidx = NULL;

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

DSP_RTN_CODE DdMWAsync::runWorkerExt()
{
	BGN_TRY_CATCH

	int signal = DSP_STAT_MW_CONTINUE;

	while (signal == DSP_STAT_MW_CONTINUE)
	{
		for (unsigned i = 0; i < worker_.size(); ++i)
		{
			switch(worker_[i]->getType())
			{
			case DdWorker::UB:
				signal = runWorkerUB(dynamic_cast<DdWorkerUB*>(worker_[i]));
				break;
			default:
				break;
			}
		}
	}

	DSPdebugMessage("runWorkerExt() terminated.\n");

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMWAsync::runWorkerUB(DdWorkerUB * workerub)
{
#define FREE_MEMORY         \
	FREE_ARRAY_PTR(number_of_elements) \
	FREE_ARRAY_PTR(indices)            \
	FREE_ARRAY_PTR(elements)

	int * number_of_elements = NULL; /**< number of elements per solution */
	int * indices            = NULL; /**< indices of solution vectors */
	double * elements        = NULL; /**< elements of solution vectors */

	/** return signal */
	int signal = DSP_STAT_MW_CONTINUE;

	BGN_TRY_CATCH

	/** upper bounds */
	vector<double> upperbounds;
	upperbounds.reserve(comm_size_);

	Solutions solutions, emptySolutions;

	DSPdebugMessage("rank %d comm_key %d\n", comm_rank_, comm_key_);
	MPI_Status status;
	int number_of_solutions;
	int recv_message = 1;

	while (recv_message)
	{
		number_of_solutions = 0;
		if (comm_key_ == 0)
		{
			MPI_Recv(&number_of_solutions, 1, MPI_INT, 0, MPI_ANY_TAG, comm_, &status);
			signal = status.MPI_TAG;
			DSPdebugMessage("number_of_solutions %d signal %d\n", number_of_solutions, signal);
		}

		/** terminate signal */
		MPI_Bcast(&signal, 1, MPI_INT, 0, subcomm_);
		if (signal == DSP_STAT_MW_STOP)
		{
			DSPdebugMessage("rank %d terminated with signal %d\n", comm_rank_, signal);
			for (unsigned i = 0; i < solutions.size(); ++i)
				FREE_PTR(solutions[i]);
			return signal;
		}

		MPI_Bcast(&number_of_solutions, 1, MPI_INT, 0, subcomm_);
		number_of_elements = new int [number_of_solutions];

		if (comm_key_ == 0)
			MPI_Recv(number_of_elements, number_of_solutions, MPI_INT, 0, MPI_ANY_TAG, comm_, &status);

		MPI_Bcast(number_of_elements, number_of_solutions, MPI_INT, 0, subcomm_);

		int total_elements = 0;
		for (int i = 0; i < number_of_solutions; ++i)
			total_elements += number_of_elements[i];
		DSPdebugMessage("total_elements %d\n", total_elements);
		indices = new int [total_elements];
		elements = new double [total_elements];

		if (comm_key_ == 0)
		{
			MPI_Recv(indices, total_elements, MPI_INT, 0, MPI_ANY_TAG, comm_, &status);
			MPI_Recv(elements, total_elements, MPI_DOUBLE, 0, MPI_ANY_TAG, comm_, &status);
			for (int i = 0, pos = 0; i < number_of_solutions; ++i)
			{
//				printf("solution %d:\n", i);
//				for (int j = 0; j < number_of_elements[i]; ++j)
//				{
//					if (j > 0 && j % 5 == 0) printf("\n");
//					printf("  [%6d] %+e", indices[pos+j], elements[pos+j]);
//				}
//				printf("\n");
				pos += number_of_elements[i];
			}
		}

		MPI_Bcast(indices, total_elements, MPI_INT, 0, subcomm_);
		MPI_Bcast(elements, total_elements, MPI_DOUBLE, 0, subcomm_);

		for (int i = 0, j = 0; i < number_of_solutions; ++i)
		{
			solutions.push_back(new CoinPackedVector(number_of_elements[i], indices + j, elements + j));
			j += number_of_elements[i];
		}

		if (comm_key_ == 0)
		{
			MPI_Iprobe(0, MPI_ANY_TAG, comm_, &recv_message, &status);
			DSPdebugMessage("MPI_Iprobe found %d.\n", recv_message);
		}
		MPI_Bcast(&recv_message, 1, MPI_INT, 0, subcomm_);
	}

	if (signal == DSP_STAT_MW_STOP)
	{
		DSPdebugMessage("rank %d terminated with signal %d\n", comm_rank_, signal);
		for (unsigned i = 0; i < solutions.size(); ++i)
			FREE_PTR(solutions[i]);
		return signal;
	}

	DSPdebugMessage("comm_key_ %d, solutions.size() %lu\n", comm_key_, solutions.size());
	if (solutions.size() > 0)
	{
		upperbounds.clear();
		for (unsigned i = 0; i < solutions.size(); ++i)
		{
			/** fix coupling solutions */
			workerub->fixCouplingVariableValues(solutions[i]);

			/** solve */
			workerub->solve();

			/** take minimum objective */
			double sumprimobj = 0.0;
			for (unsigned s = 0; s < workerub->subprobs_.size(); ++s)
				sumprimobj += workerub->subprobs_[s]->getPrimalBound();
			DSPdebugMessage("sumprimobj %e\n", sumprimobj);
			upperbounds.push_back(sumprimobj);
		}
		DSPdebugMessage("comm_key_ %d, upperbounds.size() %lu\n", comm_key_, upperbounds.size());
		double * sumub = NULL;
		if (comm_key_ == 0)
			sumub = new double [upperbounds.size()];
		MPI_Reduce(&upperbounds[0], sumub, upperbounds.size(), MPI_DOUBLE, MPI_SUM, 0, subcomm_);
		if (comm_key_ == 0)
		{
			MPI_Send(sumub, upperbounds.size(), MPI_DOUBLE, 0, 0, comm_);
			FREE_ARRAY_PTR(sumub);
		}
	}

	for (unsigned i = 0; i < solutions.size(); ++i)
		FREE_PTR(solutions[i]);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return signal;
}

DSP_RTN_CODE DdMWAsync::sendCouplingSolutions(
		int & maxnum_of_solutions,
		int & maxnum_of_indices,
		int * number_of_elements,
		int * indices,
		double * elements,
		MPI_Request * request)
{
	MPI_Status status;
	Solutions solutions; /**< solutions to distribute */

	BGN_TRY_CATCH

	DdMasterAtr * master  = dynamic_cast<DdMasterAtr*>(master_);

	/** store solutions to distribute */
	for (unsigned i = 0; i < master->worker_.size(); ++i)
	{
		for (int s = 0; s < master->nsubprobs_[i]; ++s)
		{
			int nx = model_->getNumSubproblemCouplingCols(master->subindex_[i][s]);
			CoinPackedVector * x = duplicateSolution(
					nx, master->subsolution_[i][s], ubSolutions_);
			if (x != NULL)
			{
//				DSPdebugMessage("Coupling solution:\n");
//				DSPdebug(StoMessage::printArray(nx, master->subsolution_[i][s]));
				/** store solution */
				ubSolutions_.push_back(x);
				solutions.push_back(x);
		}
	}
	}

	/** send to root process for upper-bounding */
	int number_of_solutions = solutions.size();

	if (number_of_solutions > 0)
	{
		DSPdebugMessage("send message to %d\n", comm_root_key_);
		MPI_Isend(&number_of_solutions, 1, MPI_INT, comm_root_key_, DSP_STAT_MW_CONTINUE, comm_, request);

		/** realloc if needed */
		if (maxnum_of_solutions < number_of_solutions)
		{
			number_of_elements = (int*) realloc(number_of_elements, number_of_solutions * sizeof(int));
			maxnum_of_solutions = number_of_solutions;
		}

		int total_elements = 0;
		for (int i = 0; i < number_of_solutions; ++i)
		{
			number_of_elements[i] = solutions[i]->getNumElements();
			total_elements += number_of_elements[i];
		}
		MPI_Isend(number_of_elements, number_of_solutions, MPI_INT, comm_root_key_, 0, comm_, request);

		/** realloc if needed */
		if (maxnum_of_indices < total_elements)
		{
			indices = (int*) realloc(indices, total_elements * sizeof(int));
			elements = (double*) realloc(elements, total_elements * sizeof(double));
			maxnum_of_indices = total_elements;
		}

		for (int i = 0, j = 0; i < number_of_solutions; ++i)
		{
			CoinCopyN(solutions[i]->getIndices(), solutions[i]->getNumElements(), indices + j);
			CoinCopyN(solutions[i]->getElements(), solutions[i]->getNumElements(), elements + j);
			j += solutions[i]->getNumElements();
		}
		MPI_Isend(indices, total_elements, MPI_INT, comm_root_key_, 0, comm_, request);
		MPI_Isend(elements, total_elements, MPI_DOUBLE, comm_root_key_, 0, comm_, request);
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}
