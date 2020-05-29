/*
 * DdMWAsync.cpp
 *
 *  Created on: Mar 23, 2016
 *      Author: kibaekkim
 */

//#define DSP_DEBUG
//#define DSP_DEBUG_QUEUE

#define SOLUTION_KEY_TO_STOP -999

#include "Solver/DualDecomp/DdMasterAtr.h"
#include "Solver/DualDecomp/DdMWAsync.h"
#include "Model/TssModel.h"

DdMWAsync::DdMWAsync(
		MPI_Comm     comm,   /**< MPI communicator */
		DecModel *   model,  /**< model pointer */
		DspParams *  par,    /**< parameters */
		DspMessage * message /**< message pointer */):
DdMWPara(comm,model,par,message),
qid_counter_(0),
max_queue_size_(5) {}

DdMWAsync::DdMWAsync(const DdMWAsync& rhs) :
DdMWPara(rhs),
qid_counter_(rhs.qid_counter_),
max_queue_size_(rhs.max_queue_size_),
q_id_(rhs.q_id_),
q_objval_(rhs.q_objval_) {
	double * sol;
	int * ind;
	for (int i = 0; i < rhs.q_solution_.size(); ++i) {
		// copy q solution 
		sol = new double [master_->getSiPtr()->getNumCols()];
		CoinCopyN(rhs.q_solution_[i], master_->getSiPtr()->getNumCols(), sol);
		q_solution_.push_back(sol);
		sol = NULL;

		// copy q indicator
		ind = new int [model_->getNumSubproblems()];
		CoinCopyN(rhs.q_indicator_[i], model_->getNumSubproblems(), ind);
		q_indicator_.push_back(ind);
		ind = NULL;
	}
}

DdMWAsync::~DdMWAsync()
{
	/** free queue */
	while (q_solution_.size() > 0)
		popSolutionFromQueue();
}

DSP_RTN_CODE DdMWAsync::init()
{
	BGN_TRY_CATCH

	/** This should be before DdMwPara::init(); */
	sync_ = false;

	/** initialize MPI communication settings */
	DdMWPara::init();

	if (comm_rank_ == 0)
	{
		/** create master */
		master_ = new DdMasterAtr(model_, par_, message_, lb_comm_size_);
		/** initialize master */
		master_->init();
	}
	else
	{
		/** create worker */
		if (lb_comm_rank_ >= 0)
		{
			DSPdebugMessage("Rank %d creates a worker for lower bounds.\n", comm_rank_);
			worker_.push_back(new DdWorkerLB(model_, par_, message_));
		}
		/** create worker for upper bounds */
		if (cgub_comm_rank_ >= 0)
		{
#ifdef DSP_HAS_SCIP
			DSPdebugMessage("Rank %d creates a worker for Benders cut generation.\n", comm_rank_);
			worker_.push_back(new DdWorkerCGBd(model_, par_, message_));
#endif
			DSPdebugMessage("Rank %d creates a worker for upper bounds.\n", comm_rank_);
			worker_.push_back(new DdWorkerUB(model_, par_, message_));
		}
		/** initialize workers */
		for (unsigned i = 0; i < worker_.size(); ++i)
			DSP_RTN_CHECK_THROW(worker_[i]->init());
	}
	max_queue_size_ = par_->getIntParam("DD/MAX_QSIZE");

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMWAsync::finalize()
{
	BGN_TRY_CATCH

#if 0
	char filename[64];
	const char * output_prefix = par_->getStrParam("OUTPUT/PREFIX").c_str();

	/** free master */
	if (master_)
	{
		sprintf(filename, "%s%d-DD.out", output_prefix, comm_rank_);
		writeIterInfo(filename);
		sprintf(filename, "%s%d-Master.out", output_prefix, comm_rank_);
		master_->write(filename);
		master_->finalize();
		FREE_PTR(master_);
	}

	/** free workers */
	for (unsigned i = 0; i < worker_.size(); ++i)
	{
		switch (worker_[i]->getType())
		{
		case DdWorker::LB:
			sprintf(filename, "%s%d-LB.out", output_prefix, comm_rank_);
			break;
		case DdWorker::UB:
			sprintf(filename, "%s%d-UB.out", output_prefix, comm_rank_);
			break;
		case DdWorker::CGBd:
			sprintf(filename, "%s%d-CG.out", output_prefix, comm_rank_);
			break;
		}
		worker_[i]->write(filename);
		worker_[i]->finalize();
		FREE_PTR(worker_[i]);
	}
	worker_.clear();
#else
	/** free master */
	if (master_)
	{
		master_->finalize();
		FREE_PTR(master_);
	}

	/** free workers */
	for (unsigned i = 0; i < worker_.size(); ++i)
	{
		worker_[i]->finalize();
		FREE_PTR(worker_[i]);
	}
	worker_.clear();
#endif

	/** synchronize here */
	MPI_Barrier(comm_);

	/** finalize MPI settings */
	DdMWPara::finalize();

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMWAsync::runMaster()
{
	if (comm_rank_ != 0)
		return DSP_RTN_OK;

	int signal;
	MPI_Status status;

	BGN_TRY_CATCH

	/** send start time */
	iterstime_ = CoinGetTimeOfDay();
	MPI_Bcast(&iterstime_, 1, MPI_DOUBLE, 0, comm_);

	DSP_RTN_CHECK_RTN_CODE(runMasterInit());
	DSP_RTN_CHECK_RTN_CODE(runMasterCore());

	/** TODO: Warp this part into a function, say runMasterPost(); */
	int request, dummy;
	while (1)
	{
		MPI_Iprobe(lb_comm_root_, DSP_MPI_TAG_SIG, comm_, &request, &status);
		if (request)
		{
			MPI_Recv(&signal, 1, MPI_INT, lb_comm_root_, DSP_MPI_TAG_SIG, comm_, &status);
			if (signal == DSP_STAT_MW_STOP) break;
		}
		/** CGUB worker may ask coupling solutions */
		MPI_Iprobe(cgub_comm_root_, DSP_MPI_TAG_ASK_SOLS, comm_, &request, &status);
		if (request)
		{
			/** receive dummy message */
			MPI_Recv(&dummy, 1, MPI_INT, cgub_comm_root_, DSP_MPI_TAG_ASK_SOLS, comm_, &status);
			/** send coupling solutions */
			MPIsendCoinPackedVectors(comm_, cgub_comm_root_, Solutions(), DSP_MPI_TAG_SOLS);
		}
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMWAsync::runWorker()
{
	if (comm_rank_ == 0)
		return DSP_RTN_OK;

	BGN_TRY_CATCH

	int signal = DSP_STAT_MW_CONTINUE;

	/** receive start time */
	MPI_Bcast(&iterstime_, 1, MPI_DOUBLE, 0, comm_);

	if (lb_comm_ != MPI_COMM_NULL)
	{
		DSPdebugMessage("Rank %d runs runWorkerInit() and runWorkerCore().\n", comm_rank_);
		runWorkerInit();
		runWorkerCore();

		/** TODO: Warp this part as a function, say runWorkerPost(); */
		/** make sure every processor finished */
		MPI_Barrier(lb_comm_);
		if (lb_comm_rank_ == 0)
		{
			signal = DSP_STAT_MW_STOP;
			MPI_Send(&signal, 1, MPI_INT, 0, DSP_MPI_TAG_SIG, comm_);
		}
	}
	else if (cgub_comm_ != MPI_COMM_NULL)
	{
		/** TODO: Warp this part as a function, say runWorkerCgub(); */

		int dummy = 0;
		int recv_message;
		MPI_Status status;
		Solutions solutions, local_solutions;
		while (1)
		{
			/** any signal message from lb_comm_root_ */
			if (cgub_comm_rank_ == 0)
			{
				MPI_Iprobe(lb_comm_root_, DSP_MPI_TAG_SIG, comm_, &recv_message, &status);
				DSPdebugMessage("Rank %d probed rank %d for signal (%d).\n", comm_rank_, lb_comm_root_, recv_message);
				while (recv_message)
				{
					MPI_Recv(&signal, 1, MPI_INT, lb_comm_root_, DSP_MPI_TAG_SIG, comm_, &status);
					MPI_Iprobe(lb_comm_root_, DSP_MPI_TAG_SIG, comm_, &recv_message, &status);
					DSPdebugMessage("Rank %d probed rank %d for signal (%d).\n", comm_rank_, lb_comm_root_, recv_message);
				}
			}
			/** broadcast signal */
			MPI_Bcast(&signal, 1, MPI_INT, 0, cgub_comm_);
			if (signal == DSP_STAT_MW_STOP) break;

			if (cgub_comm_rank_ == 0)
			{
				/** ask the root if it has coupling solutions to send */
				MPI_Send(&dummy, 1, MPI_INT, 0, DSP_MPI_TAG_ASK_SOLS, comm_);
				/** receive coupling solutions */
				MPIrecvCoinPackedVectors(comm_, 0, solutions, DSP_MPI_TAG_SOLS);
				message_->print(5, "The CGUB workers received %lu solutions from the root.\n", solutions.size());
			}
			/** broadcast coupling solutions */
			MPIbcastCoinPackedVectors(cgub_comm_, solutions);
			if (solutions.size() == 0) continue;

			/** TODO move some solutions to the local vector */
			for (int i = 0; i < par_->getIntParam("DD/MAX_EVAL_UB"); i++) {
				local_solutions.push_back(solutions.back());
				solutions.pop_back();
				if (solutions.size() == 0)
					break;
			}

			/** run CG */
			DSPdebugMessage("Rank %d runs runWorkerCg().\n", comm_rank_);
			signal = runWorkerCg(local_solutions);
			if (signal == DSP_STAT_MW_STOP) break;
			if (signal == DSP_STAT_MW_RESOLVE) continue;

			/** run UB */
			DSPdebugMessage("Rank %d runs runWorkerUb().\n", comm_rank_);
			signal = runWorkerUb(local_solutions);
			if (signal == DSP_STAT_MW_STOP) break;

			/** clear solutions */
			for (unsigned i = 0; i < local_solutions.size(); ++i)
				FREE_PTR(local_solutions[i]);
			local_solutions.clear();
			for (unsigned i = 0; i < solutions.size(); ++i)
				FREE_PTR(solutions[i]);
			solutions.clear();
		}
		/** clear solutions */
		for (unsigned i = 0; i < local_solutions.size(); ++i)
			FREE_PTR(local_solutions[i]);
		local_solutions.clear();
		for (unsigned i = 0; i < solutions.size(); ++i)
			FREE_PTR(solutions[i]);
		solutions.clear();

		/** send the signal receipt */
		MPI_Barrier(cgub_comm_);
		MPI_Send(&signal, 1, MPI_INT, lb_comm_root_, DSP_MPI_TAG_SIG, comm_);
		if (signal == DSP_STAT_MW_STOP && cgub_comm_rank_ == 0)
			message_->print(0, "CGUB processors are finished.\n");
	}

	DSPdebugMessage("Rank %d finished runWorker().\n", comm_rank_);

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
	}

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

	/** to store messages received from workers */
	int ** subindex = NULL;
	double ** subprimobj = NULL;
	double ** subdualobj = NULL;
	double *** subsolution = NULL;

	/** to send solutions to workers */
	Solutions solutions;

	int dummy;

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

	/** receive message */
	MPI_Gatherv(NULL, 0, MPI_DOUBLE, recvbuf, rcounts, rdispls, MPI_DOUBLE, 0, subcomm_);
#ifdef DSP_DEBUG2
	DSPdebugMessage("master receive buffer:\n");
	for (int i = 0; i < subcomm_size_; ++i)
	{
		DSPdebugMessage("  rank %d:\n", i);
		DspMessage::printArray(rcounts[i], recvbuf + rdispls[i]);
	}
#endif

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
		master->solution_key_.push_back(0);
		master->nsubprobs_.push_back(nsubprobs_[i]);
		master->subindex_.push_back(subindex[i-1]);
		master->subprimobj_.push_back(subprimobj[i-1]);
		master->subdualobj_.push_back(subdualobj[i-1]);
		master->subsolution_.push_back(subsolution[i-1]);
	}

	/** calculate dual objective */
	if (dualobj > master->bestdualobj_)
	{
		master->bestdualobj_ = dualobj;
		itercode_ = 'D';
	}

	if (parEvalUb_ >= 0 || parFeasCuts_ >= 0 || parOptCuts_ >= 0)
	{
		/** store coupling solutions */
		storeCouplingSolutions(solutions);
		/** receive dummy message */
		MPI_Recv(&dummy, 1, MPI_INT, cgub_comm_root_, DSP_MPI_TAG_ASK_SOLS, comm_, MPI_STATUS_IGNORE);
		/** send coupling solutions */
		DSP_RTN_CHECK_RTN_CODE(MPIsendCoinPackedVectors(comm_, cgub_comm_root_, solutions, DSP_MPI_TAG_SOLS));
		DSPdebugMessage("Rank %d sent %lu solutions to rank %d.\n", comm_rank_, solutions.size(), cgub_comm_root_);
		/** clear solutions */
		solutions.clear();
	}

	/** update problem */
	master->addCuts();

	/** solve problem */
	DSP_RTN_CHECK_RTN_CODE(master->solve());

	/** retrieve master solution by part */
	double * master_primsol = const_cast<double*>(master->getPrimalSolution());
	thetas  = master_primsol;
	for (int i = 0, j = model_->getNumSubproblems(); i < model_->getNumSubproblems(); ++i)
	{
		/** shallow copy */
		lambdas[i] = master_primsol + j;
		j += model_->getNumSubproblemCouplingRows(i);
	}
	/** push lambda to queue */
	if (max_queue_size_ > 0)
		DSP_RTN_CHECK_RTN_CODE(pushSolutionToQueue(master_primsol));
	for (int i = 0; i < lb_comm_size_; ++i)
		DSP_RTN_CHECK_RTN_CODE(master->setPrimsolToWorker(i, master_primsol));
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

#ifdef DSP_DEBUG2
	DSPdebugMessage("master send buffer:\n");
	for (int i = 0; i < subcomm_size_; ++i)
	{
		DSPdebugMessage("  rank %d:\n", i);
		DspMessage::printArray(scounts[i], sendbuf + sdispls[i]);
	}
#endif

	/** send message */
	MPI_Scatterv(sendbuf, scounts, sdispls, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, 0, subcomm_);

	/** release shallow-copy of pointers */
	for (int i = 0; i < model_->getNumSubproblems(); ++i)
		lambdas[i] = NULL;

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

DSP_RTN_CODE DdMWAsync::runMasterCore()
{
#define FREE_MEMORY \
	FREE_ARRAY_PTR(recvbuf) \
	FREE_2D_ARRAY_PTR(subcomm_size_-1,subindex)   \
	FREE_2D_ARRAY_PTR(subcomm_size_-1,subprimobj) \
	FREE_2D_ARRAY_PTR(subcomm_size_-1,subdualobj) \
	if (subsolution) {                                                     \
		for (int i = 0; i < subcomm_size_ - 1; ++i) {                      \
			FREE_2D_ARRAY_PTR(model_->getNumSubproblems(), subsolution[i]) \
		}                                                               \
		delete [] subsolution;                                          \
		subsolution = NULL;                                             \
	}                                                                   \
	FREE_ARRAY_PTR(numCutsAdded)

	/** MPI_Recv message:
	 *   [for each subproblem]
	 *   1 subproblem index
	 *   2 primal objective
	 *   3 dual objective
	 *   4 coupling column part of the solution
	 */
	double * recvbuf = NULL; /**< MPI_Recv: receive buffer */
	int      rcount  = 0;    /**< MPI_Recv: receive buffer size */
	int solution_key = -1;
	int nsubprobs = 0;       /**< number of subproblems to receive */

	/** to store messages received from LB workers */
	int ** subindex = NULL;
	double ** subprimobj = NULL;
	double ** subdualobj = NULL;
	double *** subsolution = NULL;

	/** MPI_Iprobe */
	int recv_message; /**< indicate if there exists a message to receive */

	int * numCutsAdded = NULL; /**< number of cuts added to each LB worker */

	/** coupling solutions */
	Solutions solutions;

	MPI_Request request_ub = MPI_REQUEST_NULL;

	/** TODO allow idle LB processors? */
	bool allowIdleLbProcessors = par_->getBoolParam("DD/ALLOW_IDLE_WORKERS");
	//int minNumWorkers

	BGN_TRY_CATCH

	int signal = DSP_STAT_MW_CONTINUE;   /**< signal to communicate with workers */
	DdMasterAtr * master = dynamic_cast<DdMasterAtr*>(master_);

	/** calculate send/receive buffer sizes */
	for (int i = 0; i < subcomm_size_; ++i)
	{
		int local_rcount = 0;
		for (int j = 0; j < nsubprobs_[i]; ++j)
			local_rcount += 3 + model_->getNumSubproblemCouplingCols(subprob_indices_[subprob_displs_[i]+j]);
		rcount = local_rcount > rcount ? local_rcount : rcount;
	}

	/** allocate memory */
	recvbuf = new double [rcount];
	subindex    = new int * [subcomm_size_-1];
	subprimobj  = new double * [subcomm_size_-1];
	subdualobj  = new double * [subcomm_size_-1];
	subsolution = new double ** [subcomm_size_-1];
	for (int i = 0; i < subcomm_size_ - 1; ++i)
	{
		subindex[i]    = new int [model_->getNumSubproblems()];
		subprimobj[i]  = new double [model_->getNumSubproblems()];
		subdualobj[i]  = new double [model_->getNumSubproblems()];
		subsolution[i] = new double * [model_->getNumSubproblems()];
		for (int s = 0; s < model_->getNumSubproblems(); ++s)
			subsolution[i][s] = new double [master_->getModelPtr()->getNumSubproblemCouplingCols(s)];
	}
	if (parFeasCuts_ >= 0 || parOptCuts_ >= 0)
	{
		numCutsAdded = new int [subcomm_size_ - 1];
		CoinZeroN(numCutsAdded, subcomm_size_ - 1);
	}

	int nTerminated = 0;
	int nIdles = 0;
	int minLbWorkers = par_->getIntParam("DD/MIN_PROCS");
	int numLbWorkers = subcomm_size_-1; /**< number of LB workers */
	double timeToRecvFirstWorker = -1.0;
	double minWaitTime = par_->getDblParam("DD/ASYNC/MIN_WAIT_TIME");

	/** print display header */
	printHeaderInfo();

	/** idle processors */
	vector<int> idle_solution_key;
	vector<int> idle_worker;
	vector<int> idle_nsubprobs;
	vector<int*> idle_subindex;

	while (signal != DSP_STAT_MW_STOP)
	{
		itercode_ = ' ';
		recv_message = 1;

		/** put idle worker processors */
		if (allowIdleLbProcessors)
		{
			nIdles = 0;
			for (unsigned i = 0; i < master->solution_key_.size(); ++i)
			{
				if (master->solution_key_[i] < 0)
				{
					idle_solution_key.push_back(-1);
					idle_worker.push_back(master->worker_[i]);
					idle_nsubprobs.push_back(master->nsubprobs_[i]);
					idle_subindex.push_back(master->subindex_[i]);
					nIdles++;
				}
			}
			message_->print(3, "Number of Idle LB workers: %d\n", nIdles);
		}
		/** clear queue */
		DSP_RTN_CHECK_RTN_CODE(master->clearSubprobData());

		timeToRecvFirstWorker = -1;

		DSPdebugMessage("################# time stamp 1: %.2f\n", CoinGetTimeOfDay() - iterstime_);
		while (1)
		{
			MPI_Status status;
			/** receive iteration signal */
			MPI_Recv(&signal, 1, MPI_INT, MPI_ANY_SOURCE, DSP_MPI_TAG_LB, subcomm_, &status);

			/** get time stamp to receive the first worker process */
			if (timeToRecvFirstWorker < 0)
				timeToRecvFirstWorker = CoinGetTimeOfDay();

			/** signal to stop? */
			if (signal == DSP_STAT_MW_STOP)
				nTerminated++;

			/** retrieve message source; the next receive should be from the same source. */
			int msg_source = status.MPI_SOURCE;
			/** receive lambda ID from the LB worker */
			MPI_Recv(&solution_key, 1, MPI_INT, msg_source, DSP_MPI_TAG_LB, subcomm_, MPI_STATUS_IGNORE);
			/** receive number of subproblems */
			MPI_Recv(&nsubprobs, 1, MPI_INT, msg_source, DSP_MPI_TAG_LB, subcomm_, MPI_STATUS_IGNORE);
			/** receive solution info from the LB worker */
			MPI_Recv(recvbuf, rcount, MPI_DOUBLE, msg_source, DSP_MPI_TAG_LB, subcomm_, MPI_STATUS_IGNORE);
#ifdef DSP_DEBUG_QUEUE
			message_->print(0,"master receive buffer (%d):\n", rcount);
			message_->print(0,"  rank %d:\n", msg_source);
			DspMessage::printArray(rcount, recvbuf);
#endif

			/** apply receive message */
			double dualobj = 0.0;
			master->solution_key_.push_back(solution_key);
			master->worker_.push_back(msg_source);
			master->nsubprobs_.push_back(nsubprobs);
			for (int s = 0, pos = 0; s < nsubprobs; ++s)
			{
				subindex[msg_source-1][s] = static_cast<int>(recvbuf[pos++]);
				subprimobj[msg_source-1][s] = recvbuf[pos++];
				subdualobj[msg_source-1][s] = recvbuf[pos++];
				CoinCopyN(recvbuf + pos,
						model_->getNumSubproblemCouplingCols(subindex[msg_source-1][s]),
						subsolution[msg_source-1][s]);
				dualobj += subprimobj[msg_source-1][s];
				pos += model_->getNumSubproblemCouplingCols(subindex[msg_source-1][s]);
				DSPdebugMessage("master received from rank %d: subprob %d primobj %+e\n",
						msg_source, subindex[msg_source-1][s], subprimobj[msg_source-1][s]);
			}
			master->subindex_.push_back(subindex[msg_source-1]);
			master->subprimobj_.push_back(subprimobj[msg_source-1]);
			master->subdualobj_.push_back(subdualobj[msg_source-1]);
			master->subsolution_.push_back(subsolution[msg_source-1]);

			/** update queue */
			for (unsigned k = 0; k < q_id_.size(); ++k) {
				if (q_id_[k] != solution_key) continue;
				q_objval_[k] += dualobj;
				for (int s = 0; s < nsubprobs; ++s)
					q_indicator_[k][subindex[msg_source-1][s]] = Q_EVALUATED;
				break;
			}

			/** check if there exists a message to receive */
			MPI_Status probe_status;
			MPI_Iprobe(MPI_ANY_SOURCE, DSP_MPI_TAG_LB, subcomm_, &recv_message, &probe_status);
#ifdef DSP_DEBUG
			if (recv_message) {
				MPI_Get_count(&probe_status, MPI_DOUBLE, &local_count);
				DSPdebugMessage("Iprobe: source %d count %d recv_message %d\n", probe_status.MPI_SOURCE, local_count, recv_message);
			}
#endif

			/** The master received messages from all the LB workers? */
			if (!recv_message &&
					(master->worker_.size() + nIdles + nTerminated >= CoinMin(numLbWorkers,minLbWorkers)
					|| minWaitTime < CoinGetTimeOfDay() - timeToRecvFirstWorker))
				break;

			/** time limit */
			if (remainingTime() < 1.0)
				break;
		}
		DSPdebugMessage("################# time stamp 2: %.2f\n", CoinGetTimeOfDay() - iterstime_);
		message_->print(3, "Number of worker processes received: %lu\n", master->worker_.size());
		if (message_->logLevel_ >= 3) {
			printf("worker processes: [");
			for (unsigned i = 0; i < master->worker_.size(); ++i)
				printf("%d ", master->worker_[i]);
			printf("]\n");
		}

		if (parEvalUb_ >= 0 || parFeasCuts_ >= 0 || parOptCuts_ >= 0)
		{
			/** store coupling solutions */
			storeCouplingSolutions(solutions);
			message_->print(3, "Added %lu solutions to the pool (%lu)\n", solutions.size(), ubSolutions_.size());

			/** probe if the worker asks solutions */
			int req_sols, dummy;
			MPI_Iprobe(cgub_comm_root_, DSP_MPI_TAG_ASK_SOLS, comm_, &req_sols, MPI_STATUS_IGNORE);
			if (req_sols)
			{
				/** receive dummy message */
				MPI_Recv(&dummy, 1, MPI_INT, cgub_comm_root_, DSP_MPI_TAG_ASK_SOLS, comm_, MPI_STATUS_IGNORE);
				/** send coupling solutions */
				MPIsendCoinPackedVectors(comm_, cgub_comm_root_, solutions, DSP_MPI_TAG_SOLS);
				/** clear solutions */
				solutions.clear();
			}

			/** receive Benders cuts */
			if (parFeasCuts_ >= 0 || parOptCuts_ >= 0)
				recvBendersCuts();
			if (parEvalUb_ >= 0)
			{
				/** Make sure that MPI_Isend for upper bound is completed. */
				if (request_ub != MPI_REQUEST_NULL)
					MPI_Wait(&request_ub, MPI_STATUS_IGNORE);
				/** receive and update upper bound */
				recvUpperBounds();
			}
		}

		/**
		 * Check the solution queues:
		 * - Check each queue from the front to the back.
		 * - Remove the queue if it is evaluated by all the worker processors.
		 * - Keep the best queue solutions among the evaluated.
		 */
		bool hasEvaluatedQueue = false;
		double * primsol_from_Q = NULL;
		double dualobj = -COIN_DBL_MAX;
		while (q_indicator_.size() > 0)
		{
			int * indicator = q_indicator_.front();
			double queue_objval = q_objval_.front();
			bool removeQueue = true;
			for (int k = 0; k < model_->getNumSubproblems(); ++k) {
				if (indicator[k] < Q_EVALUATED) {
					removeQueue = false;
					break;
				}
			}
			if (removeQueue == false) break;
			message_->print(5, "The trial point (ID %d) is evaluated and removed from the queue.\n", q_id_.front());
			hasEvaluatedQueue = true;
			if (queue_objval > dualobj)
			{
				dualobj = queue_objval;
				FREE_ARRAY_PTR(primsol_from_Q);
				primsol_from_Q = q_solution_.front();
				q_solution_.front() = NULL;
			}
			DSP_RTN_CHECK_RTN_CODE(popSolutionFromQueue());
		}

		if (hasEvaluatedQueue)
		{
			/**
			 * Update the master, if the queue with the best dual objective value is chosen.
			 */
			double olddual = master->bestdualobj_;
#ifdef DSP_DEBUG_QUEUE1
			message_->print(0, "Primal solution from Q (dual obj %+e):\n", dualobj);
			message_->printArray(master->getSiPtr()->getNumCols(), primsol_from_Q);
#endif
			message_->print(5, "best dual obj %e, current dual obj %e.\n", olddual, dualobj);
			DSP_RTN_CHECK_RTN_CODE(master->updateProblem(primsol_from_Q, dualobj));
			FREE_ARRAY_PTR(primsol_from_Q);

			/**
			 * If the best dual objective value is updated,
			 * remove all the queues not assigned to any worker processor.
			 */
			if (olddual < master->bestdualobj_)
			{
				itercode_ = itercode_ == 'P' ? 'B' : 'D';
				/** clean up queues */
				while (q_indicator_.size() > 0)
				{
					int * indicator = q_indicator_.back();
					bool removeQueue = true;
					for (int k = 0; k < model_->getNumSubproblems(); ++k)
					{
						if (indicator[k] == Q_ASSIGNED)
						{
							removeQueue = false;
							break;
						}
					}
					if (removeQueue == false) break;
					DSP_RTN_CHECK_RTN_CODE(popBackSolutionFromQueue());
				}
				message_->print(3, "The queue has been cleaned up (size %lu).\n", q_indicator_.size());
			}
		}
		else
		{
			DSP_RTN_CHECK_RTN_CODE(master->updateTrustRegion(master_->getPrimalSolution()));
		}

		/** solve problem */
		DSP_RTN_CHECK_RTN_CODE(master->solve());

		/** put solution to Q */
		if (q_solution_.size() < max_queue_size_)
		{
			/** retrieve master solution */
			double * master_primsol = const_cast<double*>(master_->getPrimalSolution());
			/** queue lambda if it is a new one. */
			DSP_RTN_CHECK_RTN_CODE(pushSolutionToQueue(master_primsol));
			master_primsol = NULL;
			message_->print(5, "Added a new trial point to the queue (size %lu).\n", q_solution_.size());
		}

		/** display iteration info */
		printIterInfo();

		/** increment iteration count */
		itercnt_++;

		/** returns continue or stop signal */
		if (remainingTime() < 1.0)
		{
			signal = DSP_STAT_MW_STOP;
			message_->print(1, "Time limit (%.2f) has been reached.\n", parTimeLimit_);
		}
		else if (itercnt_ > master_->getParPtr()->getIntParam("DD/ITER_LIM"))
		{
			signal = DSP_STAT_MW_STOP;
			message_->print(1, "Iteration limit (%d) has been reached.\n", master_->getParPtr()->getIntParam("DD/ITER_LIM"));
		} else if (nTerminated > 0) {
			signal = DSP_STAT_MW_STOP;
			message_->print(1, "One or more worker processes terminated..\n");
		} else {
			signal = master->terminationTest();
		}
		/** send signal */
		for (unsigned i = 0; i < master->worker_.size(); ++i)
		{
			DSPdebugMessage("Rank %d sent signal %d to rank %d (%d)\n", comm_rank_, signal, master->worker_[i], numLbWorkers);
			MPI_Send(&signal, 1, MPI_INT, master->worker_[i], DSP_MPI_TAG_SIG, subcomm_);
			if (signal == DSP_STAT_MW_STOP)
				message_->print(1, "The master sent STOP signal to LB processor (rank %d).\n", master->worker_[i]);
		}
		if (signal == DSP_STAT_MW_STOP)
		{
			numLbWorkers -= master->worker_.size() + idle_worker.size();
			break;
		}

		for (unsigned i = 0; i < master->worker_.size(); ++i)
		{
			/**
			 * Note that we should not see the queues evaluated here.
			 * From the first in the queue, find a queue not assigned to the worker,
			 * and assign the queue to the worker.
			 */
			double * master_primsol = NULL;
			bool chosen = chooseQueueElement(
				master->solution_key_[i], master_primsol, master->nsubprobs_[i], master->subindex_[i]);

			if (chosen == true || allowIdleLbProcessors == false) {
#ifdef DSP_DEBUG_QUEUE
				printf("number of queues %lu, master->worker_ %d, master->solution_key_ %d, q_solution_ %p\n",
						q_indicator_.size(), master->worker_[i], master->solution_key_[i], master_primsol);
#endif
				/** Get the current master solution, if no queue is assigned. */
				if (chosen == false)
					master_primsol = const_cast<double*>(master_->getPrimalSolution());

				/** Send q new trial point either from the queue or from the master */
				DSP_RTN_CHECK_RTN_CODE(
						master->setPrimsolToWorker(master->worker_[i]-1, master_primsol));
				DSP_RTN_CHECK_RTN_CODE(
						sendMasterSolution(master->solution_key_[i], master_primsol,
								master->worker_[i], master->nsubprobs_[i], master->subindex_[i],
								numCutsAdded));
			}
			master_primsol = NULL;
		}

		if (allowIdleLbProcessors)
		{
			for (unsigned i = 0; i < idle_worker.size(); ++i)
			{
				double * master_primsol = NULL;
				bool chosen = chooseQueueElement(
					idle_solution_key[i], master_primsol, idle_nsubprobs[i], idle_subindex[i]);
				if (chosen)
				{
					DSP_RTN_CHECK_RTN_CODE(master->setPrimsolToWorker(idle_worker[i]-1, master_primsol));
					DSP_RTN_CHECK_RTN_CODE(sendMasterSolution(idle_solution_key[i], master_primsol,
							idle_worker[i], idle_nsubprobs[i], idle_subindex[i], numCutsAdded));
				}
				master_primsol = NULL;
			}

			/** remove active worker processors */
			for (unsigned i = 0; i < idle_solution_key.size(); ++i)
			{
				if (idle_solution_key[i] < 0)
				{
					master->solution_key_.push_back(idle_solution_key[i]);
					master->worker_.push_back(idle_worker[i]);
					master->nsubprobs_.push_back(idle_nsubprobs[i]);
					master->subindex_.push_back(idle_subindex[i]);
				}
			}
			idle_solution_key.clear();
			idle_worker.clear();
			idle_nsubprobs.clear();
			idle_subindex.clear();
		}
	}

	/** clear solutions */
	if (parEvalUb_ >= 0 || parFeasCuts_ >= 0 || parOptCuts_ >= 0)
		solutions.clear();

	/** sending stop signals to all the worker processes */
	message_->print(1, "The master is sending STOP signal to %d LB processors.\n", numLbWorkers);
	while (numLbWorkers > 0)
	{
		MPI_Status status;
		/** receive iteration signal */
		MPI_Recv(&signal, 1, MPI_INT, MPI_ANY_SOURCE, DSP_MPI_TAG_LB, subcomm_, &status);
		/** signal to stop? */
		if (signal == DSP_STAT_MW_STOP) continue;
		/** retrieve message source; the next receive should be from the same source. */
		int msg_source = status.MPI_SOURCE;
		MPI_Recv(&solution_key, 1, MPI_INT, msg_source, DSP_MPI_TAG_LB, subcomm_, MPI_STATUS_IGNORE);
		/** receive number of subproblems */
		MPI_Recv(&nsubprobs, 1, MPI_INT, msg_source, DSP_MPI_TAG_LB, subcomm_, MPI_STATUS_IGNORE);
		/** receive solution info from the LB workers */
		MPI_Recv(recvbuf, rcount, MPI_DOUBLE, msg_source, DSP_MPI_TAG_LB, subcomm_, MPI_STATUS_IGNORE);
		/** send garbage with stop signal as a tag */
		DSPdebugMessage("Rank %d sent STOP signal to rank %d (%d)\n", comm_rank_, msg_source, numLbWorkers);
		signal = DSP_STAT_MW_STOP;
		MPI_Send(&signal, 1, MPI_INT, msg_source, DSP_MPI_TAG_SIG, subcomm_);
		message_->print(1, "The master sent STOP signal to LB processor (rank %d).\n", msg_source);
		numLbWorkers--;
	}
	for (unsigned k = 0; k < idle_worker.size(); ++k)
	{
		int dummy_signal = SOLUTION_KEY_TO_STOP;
		MPI_Send(&dummy_signal, 1, MPI_INT, idle_worker[k], DSP_MPI_TAG_LB, subcomm_);
		message_->print(0, "The master sent STOP signal to idle LB processor (rank %d).\n", idle_worker[k]);
	}
	message_->print(0, "The master is finished.\n");

	/** clear queue */
	master->clearSubprobData();

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

bool DdMWAsync::chooseQueueElement(int& qid, double*& qsol, int& nsubprobs, int*& subindex) {

	/** initialize queue ID */
	qid = -1;

	/** FIFO vs LIFO */
	for (unsigned k = 0; k < q_indicator_.size(); ++k) {
		int kk = k;
		if (par_->getBoolParam("DD/ASYNC/FIFO") == false)
			kk = q_indicator_.size() - 1 - k;
		if (q_indicator_[k][subindex[0]] == Q_NOT_ASSIGNED) {
			qid = q_id_[k];
			qsol = q_solution_[k];
			for (int s = 0; s < nsubprobs; ++s)
				q_indicator_[k][subindex[s]] = Q_ASSIGNED;
		}
		if (qid > -1) break;
	}
	DSPdebugMessage("qid %d qsol %p\n", qid, qsol);

	return (qid > -1);
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

	TssModel* tss = NULL;

	BGN_TRY_CATCH

	if (model_->isStochastic()) {
		try {
			tss = dynamic_cast<TssModel*>(model_);
		} catch (const std::bad_cast &e) {
			printf("Error: Model claims to be stochastic when it is not\n");
            return DSP_RTN_ERR;
		}
	}

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

	/** set time limit */
	workerlb->setTimeLimit(remainingTime());
	/** Solve subproblems assigned to each process  */
	workerlb->solve();

	/** create send buffer */
	for (int s = 0, pos = 0; s < narrprocidx; ++s)
	{
		sendbuf[pos++] = static_cast<double>(workerlb->subprobs_[s]->sind_);
		sendbuf[pos++] = workerlb->subprobs_[s]->getPrimalObjective();
		sendbuf[pos++] = workerlb->subprobs_[s]->getDualObjective();
		CoinCopyN(workerlb->subprobs_[s]->getSiPtr()->getColSolution(), workerlb->subprobs_[s]->ncols_coupling_, sendbuf + pos);
		pos += model_->getNumSubproblemCouplingCols(workerlb->subprobs_[s]->sind_);
		message_->print(5, "MW -> worker %d, subprob %d primobj %+e dualobj %+e\n",
				comm_rank_, workerlb->subprobs_[s]->sind_, workerlb->subprobs_[s]->getPrimalObjective(), workerlb->subprobs_[s]->getDualObjective());
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
	workerlb->solution_key_ = 0;
	for (int s = 0, pos = 0; s < narrprocidx; ++s)
	{
		workerlb->subprobs_[s]->theta_ = recvbuf[pos++];
		workerlb->subprobs_[s]->updateProblem(recvbuf + pos, tss->getProbability()[workerlb->subprobs_[s]->sind_]);
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

	/** MPI_Send message:
	 *   [for each subproblem]
	 *   1 subproblem index
	 *   2 primal objective
	 *   3 dual objective
	 *   4 coupling column part of the solution
	 */
	std::vector<double> sendbuf; /**< MPI_Send: send buffer */
	int scount = 0;              /**< MPI_Send: send buffer size */

	/** MPI_Recv message:
	 *   [for each subproblem]
	 *   1 theta
	 *   2 lambda
	 */
	std::vector<double> recvbuf; /**< MPI_Recv: receive buffer */
	int rcount = 0;              /**< MPI_Recv: receive buffer size */

	/** separate MPI_Recv */
	int nsubprobs = 0;         /**< number of subproblems to receive */
	std::vector<int> subindex; /**< subproblem indices to receive */

	BGN_TRY_CATCH

	assert(worker_[0]->getType()==DdWorker::LB);
	DdWorkerLB * workerlb = dynamic_cast<DdWorkerLB*>(worker_[0]);

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
	sendbuf.resize(scount);
	recvbuf.resize(rcount);

	/** loop until when the master signals stop */
	while(1)
	{
#if 0
		/** set gap tolerance */
		double gaptol = 100.0;
		if (signal == DSP_STAT_MW_EXACT)
			gaptol = par_->getDblParam("DD/STOP_TOL");
		for (int i = 0; i < narrprocidx; ++i)
			workerlb->subprobs_[i]->setGapTol(gaptol);
#endif

		/** set time limit */
		workerlb->setTimeLimit(remainingTime());
		/** Solve subproblems assigned to each process  */
		workerlb->solve();

		/** worker status */
		signal = workerlb->getStatus();
		/** send stop signal to master */
		MPI_Send(&signal, 1, MPI_INT, 0, DSP_MPI_TAG_LB, subcomm_);

		if (signal != DSP_STAT_MW_STOP) {
			/** create send buffer */
			for (int s = 0, pos = 0; s < workerlb->getNumSubprobs(); ++s)
			{
				sendbuf[pos++] = static_cast<double>(workerlb->subprobs_[s]->sind_);
				sendbuf[pos++] = workerlb->subprobs_[s]->getPrimalObjective();
				sendbuf[pos++] = workerlb->subprobs_[s]->getDualObjective();
				CoinCopyN(workerlb->subprobs_[s]->getSiPtr()->getColSolution(), workerlb->subprobs_[s]->ncols_coupling_, &sendbuf[0] + pos);
				pos += workerlb->subprobs_[s]->ncols_coupling_;
				DSPdebugMessage("worker %d, subprob %d primobj %+e dualobj %+e\n",
						comm_rank_, workerlb->subprobs_[s]->sind_, workerlb->subprobs_[s]->getPrimalObjective(), workerlb->subprobs_[s]->getDualObjective());
			}
#ifdef DSP_DEBUG_QUEUE1
			message_->print(0, "LB processor (rank %d) send message (%d):\n", comm_rank_, scount);
			DspMessage::printArray(scount, &sendbuf[0]);
#endif

			/** send lambda ID to the master */
			MPI_Send(&(workerlb->solution_key_), 1, MPI_INT, 0, DSP_MPI_TAG_LB, subcomm_);
//			message_->print(1, "LB processor (rank %d) evaluated the trial point (ID %d).\n", comm_rank_, workerlb->solution_key_);
			/** send number of subproblems to the master */
			nsubprobs = workerlb->getNumSubprobs();
			MPI_Send(&nsubprobs, 1, MPI_INT, 0, DSP_MPI_TAG_LB, subcomm_);
			/** send message to the master */
			MPI_Send(&sendbuf[0], scount, MPI_DOUBLE, 0, DSP_MPI_TAG_LB, subcomm_);
		}
		/** receive signal from the root */
		MPI_Recv(&signal, 1, MPI_INT, 0, DSP_MPI_TAG_SIG, subcomm_, MPI_STATUS_IGNORE);
		if (signal == DSP_STAT_MW_STOP)
		{
			DSPdebugMessage("LB processor (rank %d) received STOP signal.\n", comm_rank_);
			break;
		}

		/** receive lambda ID from the master */
		MPI_Recv(&(workerlb->solution_key_), 1, MPI_INT, 0, DSP_MPI_TAG_LB, subcomm_, MPI_STATUS_IGNORE);
		if (workerlb->solution_key_ == SOLUTION_KEY_TO_STOP) break;
		//message_->print(1, "LB processor (rank %d) received a trial point (ID %d).\n", comm_rank_, workerlb->solution_key_);

		/** receive number of subproblms */
		MPI_Recv(&nsubprobs, 1, MPI_INT, 0, DSP_MPI_TAG_LB, subcomm_, MPI_STATUS_IGNORE);

		/** receive subproblm indices */
		subindex.resize(nsubprobs);
		MPI_Recv(&subindex[0], nsubprobs, MPI_INT, 0, DSP_MPI_TAG_LB, subcomm_, MPI_STATUS_IGNORE);

		/** calculate size of send/receive buffer */
		scount = 0; rcount = 0;
		for (int i = 0; i < nsubprobs; ++i) {
			scount += 3 + model_->getNumSubproblemCouplingRows(subindex[i]);
			rcount += 1 + model_->getNumSubproblemCouplingRows(subindex[i]);
		}
		sendbuf.resize(scount);
		recvbuf.resize(rcount);

		/** receive message from the master */
		MPI_Recv(&recvbuf[0], rcount, MPI_DOUBLE, 0, DSP_MPI_TAG_LB, subcomm_, MPI_STATUS_IGNORE);
#ifdef DSP_DEBUG_QUEUE1
		if (comm_rank_ == 1) {
			printf("LB worker (rank %d) received message (%d):\n", comm_rank_, rcount);
			DspMessage::printArray(rcount, &recvbuf[0]);
		}
#endif
		/** receive cuts */
		if (parFeasCuts_ >= 0 || parOptCuts_ >= 0) {
			/** prepare cuts to receive */
			OsiCuts cutsToRecv;
			/** receive cuts from the master */
			MPIrecvOsiCuts(subcomm_, 0, cutsToRecv, DSP_MPI_TAG_CGBD);
			DSPdebugMessage("Rank %d received %d cuts from rank 0.\n", comm_rank_, cutsToRecv.sizeCuts());
			/** move cuts */
			for (int i = 0; i < cutsToRecv.sizeCuts(); ++i) {
				OsiRowCut * rc = cutsToRecv.rowCutPtr(i);
				cutsToAdd_->insert(rc);
			}
			cutsToRecv.dumpCuts();
		}

		/** receive upper bounds */
		double bestprimobj = COIN_DBL_MAX;
		if (parEvalUb_ >= 0) {
			MPI_Recv(&bestprimobj, 1, MPI_DOUBLE, 0, DSP_MPI_TAG_UB, subcomm_, MPI_STATUS_IGNORE);
			DSPdebugMessage("Rank %d received upper bound %+e from rank 0.\n", comm_rank_, bestprimobj);
		}

		/** set workerlb */
		DSP_RTN_CHECK_THROW(setWorkerLb(workerlb, nsubprobs, &subindex[0], &recvbuf[0], bestprimobj));
	}

	if (lb_comm_rank_ == 0) {
		message_->print(0, "LB processors are finished.\n");
		if (parFeasCuts_ >= 0 || parOptCuts_ >= 0 || parEvalUb_ >= 0) {
			message_->print(0, "Finishing CGUB processors...\n");
			MPI_Send(&signal, 1, MPI_INT, cgub_comm_root_, DSP_MPI_TAG_SIG, comm_);
			message_->print(3, "The LB root worker (rank %d) sent STOP signal to the CGUB root worker (rank %d).\n",
					comm_rank_, cgub_comm_root_);
			MPI_Recv(&signal, 1, MPI_INT, cgub_comm_root_, DSP_MPI_TAG_SIG, comm_, MPI_STATUS_IGNORE);
			message_->print(3, "The LB root worker (rank %d) received the STOP signal receipt from the CGUB root worker (rank %d).\n",
					comm_rank_, cgub_comm_root_);
		}
	}

	/** release pointers */
	arrprocidx = NULL;

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMWAsync::setWorkerLb(DdWorkerLB* workerlb, int nsubprobs, int* subindex, double* buf, double bestprimobj) {

	TssModel* tss = NULL;

	BGN_TRY_CATCH

	if (model_->isStochastic()) {
		try {
			tss = dynamic_cast<TssModel*>(model_);
		} catch (const std::bad_cast &e) {
			printf("Error: Model claims to be stochastic when it is not\n");
            return DSP_RTN_ERR;
		}
	}

	for (int s = 0, pos = 0; s < nsubprobs; ++s) {
		workerlb->subprobs_[s]->theta_ = buf[pos++];
		workerlb->subprobs_[s]->updateProblem(buf + pos, tss->getProbability()[workerlb->subprobs_[s]->sind_], bestprimobj);
		/** apply Benders cuts */
		if (cutsToAdd_->sizeCuts() > 0 && (parFeasCuts_ >= 0 || parOptCuts_ >= 0)) {
			workerlb->subprobs_[s]->pushCuts(cutsToAdd_);
			DSPdebugMessage("Rank %d pushed %d Benders cuts.\n", comm_rank_, cutsToAdd_->sizeCuts());
		}
		pos += model_->getNumSubproblemCouplingRows(workerlb->subprobs_[s]->sind_);
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMWAsync::storeCouplingSolutions(Solutions& solutions) {

	BGN_TRY_CATCH

	DdMasterAtr * master  = dynamic_cast<DdMasterAtr*>(master_);

	/** store solutions to distribute */
	for (unsigned i = 0; i < master->nsubprobs_.size(); ++i)
		for (int s = 0; s < master->nsubprobs_[i]; ++s)
		{
			int nx = model_->getNumSubproblemCouplingCols(master->subindex_[i][s]);
			DSPdebug2({
				DSPdebugMessage("ubSolutions_ %lu\n", ubSolutions_.size());
				DSPdebugMessage("solution[%d] nx %d:\n", s, nx);
				message_->printArray(nx, master->subsolution_[i][s]);});

			CoinPackedVector * x = duplicateSolution(
					nx, master->subsolution_[i][s], ubSolutions_);
			if (x != NULL)
			{
				DSPdebug2({
					DSPdebugMessage2("Coupling solution:\n");
					DspMessage::printArray(nx, master->subsolution_[i][s]);});
				/** store solution */
				ubSolutions_.push_back(x);
				solutions.push_back(x);
			}
		}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMWAsync::runWorkerCg(
		Solutions solutions /**< solutions at which cuts are generated */)
{
	if (parFeasCuts_ < 0 && parOptCuts_ < 0) return DSP_RTN_OK;
	DSPdebugMessage("Rank %d cgub_comm_rank_ %d\n", comm_rank_, cgub_comm_rank_);

	OsiCuts cuts;
	int signal = DSP_STAT_MW_CONTINUE;

	BGN_TRY_CATCH
#ifdef DSP_HAS_SCIP
	if (solutions.size() > 0)
	{
		/** generate Benders cuts */
		int cg_status = generateBendersCuts(cgub_comm_, cgub_comm_rank_, solutions, cuts);
		MPI_Allreduce(&cg_status, &signal, 1, MPI_INT, MPI_MAX, cgub_comm_);
		/** FIXME: This must be wrong. All cgub processors should send cuts.
		 * See the corresponding part in the master. */
		//if (cgub_comm_rank_ == 0)
		{
			DSPdebugMessage("Rank %d has cut generation status %d.\n", comm_rank_, signal);
			/** send the cuts to the root */
			MPIsendOsiCuts(comm_, 0, cuts, DSP_MPI_TAG_CGBD);
			DSPdebugMessage("Rank %d sent %d cuts to the root.\n", comm_rank_, cuts.sizeCuts());
		}
	}
#endif
	END_TRY_CATCH_RTN(;,DSP_STAT_MW_STOP)

	return signal;
}

DSP_RTN_CODE DdMWAsync::runWorkerUb(
		Solutions solutions /**< solutions to evaluate UB */) {
	if (parEvalUb_ < 0) return DSP_RTN_OK;
	DSPdebugMessage("rank %d ub_comm_rank_ %d\n", comm_rank_, cgub_comm_rank_);

	DdWorkerUB * workerub = NULL;
	vector<double> upperbounds;
	/** default signal */
	int signal = DSP_STAT_MW_CONTINUE;

	BGN_TRY_CATCH

	/** retrieve upper bounding worker */
	for (unsigned i = 0; i < worker_.size(); ++i)
	{
		if (worker_[i]->getType() == DdWorker::UB)
		{
			workerub = dynamic_cast<DdWorkerUB*>(worker_[i]);
			DSPdebugMessage("Rank %d works for upper bounds.\n", comm_rank_);
			break;
		}
	}
	if (workerub == NULL)
	{
		DSPdebugMessage("Rank %d does not work for upper bounds.\n", comm_rank_);
		return DSP_STAT_MW_STOP;
	}

	//message_->print(0, "W%d  Starting UB for %lu solutions.\n", comm_rank_, solutions.size());
	if (solutions.size() > 0)
	{
#if 0
		if (cgub_comm_rank_ == 0) {
			char filename[64];
			const char * output_prefix = par_->getStrParam("OUTPUT/PREFIX").c_str();
			sprintf(filename, "%s-solns.txt", output_prefix);

			std::ofstream outfile;
			outfile.open(filename, std::ios_base::app);
			for (unsigned i = 0; i < solutions.size(); ++i) {
				CoinPackedVector* soln = solutions[i];
				bool bewritten = false;
				for (int j = 0; j < soln->getNumElements(); ++j) {
					if (fabs(soln->getElements()[j]) < 1.0e-8) continue;
					if (bewritten) outfile << ",";
					outfile << soln->getIndices()[j];
					bewritten = true;
				}
				outfile << "\n";
			}
			outfile.close();
		}
#endif
		upperbounds.reserve(comm_size_);
		upperbounds.clear();
		for (unsigned i = 0; i < solutions.size(); ++i)
		{
			/** set time limit */
			workerub->setTimeLimit(remainingTime());
			/** evaluate upper bounds */
			double sumprimobj = workerub->evaluate(solutions[i]);
			DSPdebugMessage("Rank %d: sumprimobj %e\n", comm_rank_, sumprimobj);
			upperbounds.push_back(sumprimobj);
		}
		DSPdebugMessage("Rank %d: cgub_comm_rank_ %d, upperbounds.size() %lu\n", comm_rank_, cgub_comm_rank_, upperbounds.size());
		double * sumub = NULL;
		if (cgub_comm_rank_ == 0)
			sumub = new double [upperbounds.size()];
		MPI_Reduce(&upperbounds[0], sumub, upperbounds.size(), MPI_DOUBLE, MPI_SUM, 0, cgub_comm_);
		if (cgub_comm_rank_ == 0)
		{
			bool updated = false;
			for (unsigned i = 0; i < upperbounds.size(); ++i)
				if (workerub->bestub_ > sumub[i]) {
					workerub->bestub_ = sumub[i];
					updated = true;
				}
			if (updated)
				MPI_Send(&(workerub->bestub_), 1, MPI_DOUBLE, 0, DSP_MPI_TAG_UB, comm_);
			FREE_ARRAY_PTR(sumub);
		}
	}
	//message_->print(0, "W%d  Ended UB for %lu solutions.\n", comm_rank_, solutions.size());

	END_TRY_CATCH_RTN(;,DSP_STAT_MW_STOP)

	return signal;
}
#if 0
DSP_RTN_CODE DdMWAsync::recvCouplingSolutions(
		MPI_Comm comm, /**< communicator to broadcast solutions */
		int comm_rank, /**< processor rank of the given communicator */
		Solutions &solutions /**< received solution placeholder */) {
#define FREE_MEMORY         \
	FREE_ARRAY_PTR(number_of_elements) \
	FREE_ARRAY_PTR(indices)            \
	FREE_ARRAY_PTR(elements)

	int * number_of_elements = NULL; /**< number of elements per solution */
	int * indices            = NULL; /**< indices of solution vectors */
	double * elements        = NULL; /**< elements of solution vectors */

	MPI_Status status;
	int flag;
	int number_of_solutions; /** number of solutions to receive */

	BGN_TRY_CATCH

	/** clear the local solution pool */
	for (unsigned i = 0; i < solutions.size(); ++i)
		FREE_PTR(solutions[i]);
	solutions.clear();

	/** receive the number of solutions */
	if (comm_rank == 0)
	{
		MPI_Recv(&number_of_solutions, 1, MPI_INT, 0, DSP_MPI_TAG_SOLS, comm_, &status);
		DSPdebugMessage("Rank %d received the number_of_solutions %d\n", comm_rank_, number_of_solutions);
	}
	/** broadcast the number of solutions */
	MPI_Bcast(&number_of_solutions, 1, MPI_INT, 0, comm);

	/** receive the number of elements for each solution */
	number_of_elements = new int [number_of_solutions];
	if (comm_rank == 0)
		MPI_Recv(number_of_elements, number_of_solutions, MPI_INT, 0, DSP_MPI_TAG_SOLS, comm_, &status);
	/** broadcast the number of elements */
	MPI_Bcast(number_of_elements, number_of_solutions, MPI_INT, 0, comm);

	int total_elements = 0;
	for (int i = 0; i < number_of_solutions; ++i)
		total_elements += number_of_elements[i];
	DSPdebugMessage("Rank %d received total_elements %d\n", comm_rank_, total_elements);

	/** receive indices and elements */
	indices = new int [total_elements];
	elements = new double [total_elements];
	if (comm_rank == 0)
	{
		MPI_Recv(indices, total_elements, MPI_INT, 0, DSP_MPI_TAG_SOLS, comm_, &status);
		MPI_Recv(elements, total_elements, MPI_DOUBLE, 0, DSP_MPI_TAG_SOLS, comm_, &status);
	}
	/** broadcast indices and elements */
	MPI_Bcast(indices, total_elements, MPI_INT, 0, comm);
	MPI_Bcast(elements, total_elements, MPI_DOUBLE, 0, comm);

	for (int i = 0, j = 0; i < number_of_solutions; ++i)
	{
		solutions.push_back(new CoinPackedVector(number_of_elements[i], indices + j, elements + j));
		j += number_of_elements[i];
	}
	DSPdebugMessage("Rank %d has solutions.size() %lu\n", comm_rank_, solutions.size());

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}
#endif
DSP_RTN_CODE DdMWAsync::recvBendersCuts() {

	MPI_Status status;
	int recv_message;
	int cgub_proc;
	OsiCuts cuts;

	BGN_TRY_CATCH

	/** FIXME: The master should receive cuts from all cgub processors.
	 * See the corresponding part in the worker. */

	/** is there a message to receive? */
	MPI_Iprobe(MPI_ANY_SOURCE, DSP_MPI_TAG_CGBD, comm_, &recv_message, &status);
	//MPI_Iprobe(cgub_comm_root_, DSP_MPI_TAG_CGBD, comm_, &recv_message, &status);
	//DSPdebugMessage("MPI_Iprobe returns recv_message %d from rank %d.\n" ,recv_message, cgub_comm_root_);
	if (recv_message == 0) return DSP_RTN_ERR;

	/** receive cuts from CG worker */
	while (recv_message)
	{
		/** identify the source */
		cgub_proc = status.MPI_SOURCE;
		DSPdebugMessage("MPI_Iprobe returns recv_message %d from rank %d.\n" ,recv_message, cgub_proc);
		/** receive cuts */
		MPIrecvOsiCuts(comm_, cgub_proc, cuts, DSP_MPI_TAG_CGBD);
		//MPIrecvOsiCuts(comm_, cgub_comm_root_, cuts, DSP_MPI_TAG_CGBD);
		/** is there a message to receive? */
		MPI_Iprobe(MPI_ANY_SOURCE, DSP_MPI_TAG_CGBD, comm_, &recv_message, &status);
		//MPI_Iprobe(cgub_comm_root_, DSP_MPI_TAG_CGBD, comm_, &recv_message, &status);
	}

	/** move cuts */
	for (int i = 0; i < cuts.sizeCuts(); ++i)
	{
		OsiRowCut * rc = cuts.rowCutPtr(i);
		cutsToAdd_->insert(rc);
	}
	cuts.dumpCuts();

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMWAsync::recvUpperBounds() {
#define FREE_MEMORY \
	FREE_ARRAY_PTR(primobjs)

	MPI_Status status;
	int recv_message;
	double * primobjs = NULL; /**< primal objectives from subproblems */
	int nprimobjs; /**< number of UBs to receive */

	BGN_TRY_CATCH

	/** is there a message to receive? */
	MPI_Iprobe(cgub_comm_root_, DSP_MPI_TAG_UB, comm_, &recv_message, &status);
	DSPdebugMessage("MPI_Iprobe returns recv_message %d.\n" ,recv_message);
	if (recv_message == 0) return DSP_RTN_ERR;

	/** allocate memory */
	primobjs = new double [model_->getNumSubproblems()];

	/** receive upper bounds and cuts */
	double bestprimobj = master_->bestprimobj_;
	while (recv_message)
	{
		/** get number of UBs to receive */
		MPI_Get_count(&status, MPI_DOUBLE, &nprimobjs);
		/** receive UBs */
		MPI_Recv(primobjs, nprimobjs, MPI_DOUBLE, cgub_comm_root_, DSP_MPI_TAG_UB, comm_, &status);
		/** calculate the best UB */
		for (int i = 0; i < nprimobjs; ++i)
		{
			message_->print(5, "solution %d: primal objective %+e\n", i, primobjs[i]);
			master_->bestprimobj_ = primobjs[i] < master_->bestprimobj_ ? primobjs[i] : master_->bestprimobj_;
		}
		/** is there a message to receive? */
		MPI_Iprobe(cgub_comm_root_, DSP_MPI_TAG_UB, comm_, &recv_message, &status);
	}

	if (bestprimobj > master_->bestprimobj_)
		itercode_ = itercode_ == 'D' ? 'B' : 'P';

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

DSP_RTN_CODE DdMWAsync::sendMasterSolution(
		int solution_key,
		double * master_primsol,
		int worker_proc,
		int num_subprobs,
		int * subprobs,
		int * numCutsAdded)
{
#define FREE_MEMORY \
	FREE_ARRAY_PTR(sendbuf)

	/** MPI_Send message:
	 *   [for each subproblem]
	 *   1 theta
	 *   2 lambda
	 */
	double * sendbuf = NULL; /**< MPI_Send: send buffer */
	int      scount  = 0;    /**< MPI_Send: send buffer size */
	double * thetas  = NULL; /**< theta part of the solution */
	double * lambdas = NULL; /**< lambda part of the solution */
	int local_scount;

	BGN_TRY_CATCH

	/** calculate send/receive buffer sizes */
	for (int i = 0; i < subcomm_size_; ++i)
	{
		local_scount = 0;
		for (int j = 0; j < nsubprobs_[i]; ++j)
			local_scount += 1 + model_->getNumSubproblemCouplingRows(subprob_indices_[subprob_displs_[i]+j]);
		scount = local_scount > scount ? local_scount : scount;
	}

	/** allocate memory */
	sendbuf = new double [scount];

	/** create send buffer */
	thetas = master_primsol;
	lambdas = master_primsol + model_->getNumSubproblems();
	local_scount = 0;
	int lambda_offset = 0;
	for (unsigned s = 0, k = 0; s < num_subprobs;)
	{
		if (k < subprobs[s])
		{
			lambdas += model_->getNumSubproblemCouplingRows(k);
			lambda_offset += model_->getNumSubproblemCouplingRows(k);
			k++;
			continue;
		}
		sendbuf[local_scount++] = thetas[subprobs[s]];
		CoinCopyN(lambdas,
				model_->getNumSubproblemCouplingRows(subprobs[s]),
				sendbuf + local_scount);
		local_scount += model_->getNumSubproblemCouplingRows(subprobs[s]);
		s++;
	}

	/** send lambda ID */
	MPI_Send(&(solution_key), 1, MPI_INT, worker_proc, DSP_MPI_TAG_LB, subcomm_);

	/** send number of subproblems */
	MPI_Send(&num_subprobs, 1, MPI_INT, worker_proc, DSP_MPI_TAG_LB, subcomm_);

	/** send subproblem indices */
	MPI_Send(subprobs, num_subprobs, MPI_INT, worker_proc, DSP_MPI_TAG_LB, subcomm_);

	/** send message */
	MPI_Send(sendbuf, local_scount, MPI_DOUBLE, worker_proc, DSP_MPI_TAG_LB, subcomm_);
#ifdef DSP_DEBUG_QUEUE1
	message_->print(0, "The master sent a trial point (ID %d) to LB processor (rank %d).\n", solution_key, worker_proc);
	DspMessage::printArray(local_scount, sendbuf);
#endif
	if (parFeasCuts_ >= 0 || parOptCuts_ >= 0)
	{
		/** prepare cuts to send */
		OsiCuts cutsToSend;
		for (int j = numCutsAdded[worker_proc-1]; j < cutsToAdd_->sizeCuts(); ++j)
		{
			OsiRowCut * rc = cutsToAdd_->rowCutPtr(j);
			cutsToSend.insert(rc);
		}
		numCutsAdded[worker_proc-1] = cutsToSend.sizeCuts();
		/** send cuts to LB workers */
		MPIsendOsiCuts(subcomm_, worker_proc, cutsToSend, DSP_MPI_TAG_CGBD);
		DSPdebugMessage("Rank %d sent %d cuts to rank %d.\n", comm_rank_, cutsToSend.sizeCuts(), worker_proc);
		cutsToSend.dumpCuts();
	}

	/** send upper bounds */
	if (parEvalUb_ >= 0)
	{
		MPI_Send(&(master_->bestprimobj_), 1, MPI_DOUBLE, worker_proc, DSP_MPI_TAG_UB, subcomm_);
		DSPdebugMessage("Rank %d sent upper bound %+e to rank %d.\n", comm_rank_, master_->bestprimobj_, worker_proc);
	}

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY
	thetas = NULL;

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMWAsync::pushSolutionToQueue(double* solution)
{
	BGN_TRY_CATCH

	/** copy memory */
	double * l = new double [master_->getSiPtr()->getNumCols()];
	int * indicator = new int [model_->getNumSubproblems()];
	CoinCopyN(solution, master_->getSiPtr()->getNumCols(), l);
	CoinFillN(indicator, model_->getNumSubproblems(), 0);

	/** push data */
	q_id_.push_back(qid_counter_++);
	q_solution_.push_back(l);
	q_indicator_.push_back(indicator);
	q_objval_.push_back(0.0);
	l = NULL;
	indicator = NULL;

#ifdef DSP_DEBUG_QUEUE1
	DSPdebugMessage("Print Queue:\n");
	for (unsigned i = 0; i < q_id_.size(); ++i)
		printf("  q_id_ %d, q_solution_ %p, q_objval_ %e\n", q_id_[i], q_solution_[i], q_objval_[i]);
#endif

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMWAsync::pushFrontSolutionToQueue(double* solution)
{
	BGN_TRY_CATCH

	/** copy memory */
	double * l = new double [master_->getSiPtr()->getNumCols()];
	int * indicator = new int [model_->getNumSubproblems()];
	CoinCopyN(solution, master_->getSiPtr()->getNumCols(), l);
	CoinFillN(indicator, model_->getNumSubproblems(), 0);

	/** push data */
	q_id_.push_front(qid_counter_++);
	q_solution_.push_front(l);
	q_indicator_.push_front(indicator);
	q_objval_.push_front(0.0);
	l = NULL;
	indicator = NULL;

#ifdef DSP_DEBUG_QUEUE1
	DSPdebugMessage("Print Queue:\n");
	for (unsigned i = 0; i < q_id_.size(); ++i)
		printf("  q_id_ %d, q_solution_ %p, q_objval_ %e\n", q_id_[i], q_solution_[i], q_objval_[i]);
#endif

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMWAsync::popSolutionFromQueue()
{
	BGN_TRY_CATCH

	if (q_id_.size() <= 0) throw "Queue q_id_ is empty.";
	if (q_solution_.size() <= 0) throw "Queue q_lambda_ is empty.";
	if (q_indicator_.size() <= 0) throw "Queue q_indicator_ is empty.";
	if (q_objval_.size() <= 0) throw "Queue q_objval_ is empty.";

	/** free memeory */
	FREE_ARRAY_PTR(q_solution_.front());
	FREE_ARRAY_PTR(q_indicator_.front());

	q_id_.pop_front();
	q_solution_.pop_front();
	q_indicator_.pop_front();
	q_objval_.pop_front();

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMWAsync::popBackSolutionFromQueue()
{
	BGN_TRY_CATCH

	if (q_id_.size() <= 0) throw "Queue q_id_ is empty.";
	if (q_solution_.size() <= 0) throw "Queue q_lambda_ is empty.";
	if (q_indicator_.size() <= 0) throw "Queue q_indicator_ is empty.";
	if (q_objval_.size() <= 0) throw "Queue q_objval_ is empty.";

	/** free memeory */
	FREE_ARRAY_PTR(q_solution_.back());
	FREE_ARRAY_PTR(q_indicator_.back());

	q_id_.pop_back();
	q_solution_.pop_back();
	q_indicator_.pop_back();
	q_objval_.pop_back();

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}
