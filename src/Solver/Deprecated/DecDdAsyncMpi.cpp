/*
 * DecDdAsyncMpi.cpp
 *
 *  Created on: Feb 4, 2016
 *      Author: kibaekkim
 */

#include <Solver/DecDdAsyncMpi.h>

DecDdAsyncMpi::DecDdAsyncMpi(MPI_Comm comm) :
DecDdMpi(comm, "DecDdAsyncMip.log"),
scount_(0), rcount_(0),
sendbuf_(NULL), recvbuf_(NULL),
maxnumsubprobs_(0), subindex_(NULL), subprimobj_(NULL), subdualobj_(NULL), subsolution_(NULL)
{
	/** indicate that the master is dedicated to the master problem, not solving any subproblem. */
	use_root_process_ = false;
}

DecDdAsyncMpi::~DecDdAsyncMpi()
{
	FREE_ARRAY_PTR(sendbuf_);
	FREE_ARRAY_PTR(recvbuf_);
	FREE_ARRAY_PTR(recvbuf_);
	FREE_ARRAY_PTR(subindex_);
	FREE_ARRAY_PTR(subprimobj_);
	FREE_ARRAY_PTR(subdualobj_);
	FREE_2D_ARRAY_PTR(maxnumsubprobs_,subsolution_);
}

DSP_RTN_CODE DecDdAsyncMpi::solve()
{
	if (comm_size_ < 2)
	{
		message_->print(0, "Error: At least two processors are required.\n");
		status_ = STO_STAT_STOPPED_MPI;
		return STO_RTN_ERR;
	}

	BGN_TRY_CATCH

	/** initialize global settings */
	initializeGlobal();

	/** initialize local settings */
	initializeLocal();

	/** create master problem */
	createMaster();

	/** create subproblems */
	createSubproblem();

	/** initialize problem-dependent statistics */
	initializeStatistics();

	/** run master process */
	runMaster();

	/** run worker processes */
	runWorkers();

	/** finalize */

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
}

/** initialize local settings */
DSP_RTN_CODE DecDdAsyncMpi::initializeLocal()
{
	int nsubprobs = 0; /**< number of subproblems */

	BGN_TRY_CATCH

	/** determine the sizes of sending and receiving messages. */
	if (comm_rank_ != 0)
	{
		/** sending message size for worker
		 * number of subproblems
		 * [for each subproblem]
		 *   subproblem index
		 *   primal objective
		 *   dual objective
		 *   subproblem solution for coupling columns
		 * */
		scount_ = 1;
		for (int s = 0; s < parProcIdxSize_; ++s)
			scount_ += 3 + model_->getNumSubproblemCouplingCols(parProcIdx_[s]);
		/** receiving message size for worker
		 * signal (continue or stop?)
		 * theta
		 * lambda
		 * */
		rcount_ = 1 + parProcIdxSize_;// + model_->getNumCouplingRows();
		for (int s = 0; s < parProcIdxSize_; ++s)
			rcount_ += model_->getNumSubproblemCouplingRows(s);
		/** number of subproblems */
		nsubprobs = parProcIdxSize_;
	}

	int tmprcount = 0;
	/** receiving message size for master */
	MPI_Reduce(&scount_, &tmprcount, 1, MPI_INT, MPI_MAX, 0, comm_);
	/** sending message size for master */
	MPI_Reduce(&rcount_, &scount_, 1, MPI_INT, MPI_MAX, 0, comm_);
	if (comm_rank_ == 0)
		rcount_ = tmprcount;
	/** maximum number of subproblems */
	MPI_Reduce(&nsubprobs, &maxnumsubprobs_, 1, MPI_INT, MPI_MAX, 0, comm_);

	/** allocate memory */
	sendbuf_     = new double [scount_];
	recvbuf_     = new double [rcount_];
	if (comm_rank_ == 0)
	{
		subindex_    = new int [maxnumsubprobs_];
		subprimobj_  = new double [maxnumsubprobs_];
		subdualobj_  = new double [maxnumsubprobs_];
		subsolution_ = new double * [maxnumsubprobs_];
		for (int s = 0; s < maxnumsubprobs_; ++s)
			subsolution_[s] = new double [model_->getNumSubproblemCouplingCols(s)];
	}

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
}

/** create master problem */
DSP_RTN_CODE DecDdAsyncMpi::createMaster()
{
	if (comm_rank_ != 0)
		return STO_RTN_OK;

	BGN_TRY_CATCH

	/** TODO: set master method */

	/** create master problem */
	master_->createProblem(model_);

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
}

/** TODO: run master process */
DSP_RTN_CODE DecDdAsyncMpi::runMaster()
{
	if (comm_rank_ != 0)
		return STO_RTN_OK;

	MPI_Status status;
	int nworkers = comm_size_ - 1; /**< number of workers */
	int nsubprobs = 0; /**< number of subproblems received from worker */
	int pos = 0;
	DSP_RTN_CODE master_signal = STO_STAT_MW_CONTINUE;

	BGN_TRY_CATCH

	while (master_signal == STO_STAT_MW_CONTINUE)
	{
		/** receive message */
		MPI_Recv(recvbuf_, rcount_, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, comm_, &status);
		nsubprobs = static_cast<int>(recvbuf_[pos++]);
		for (int s = 0; s < nsubprobs; ++s)
		{
			subindex_[s] = static_cast<int>(recvbuf_[pos++]);
			subprimobj_[s] = recvbuf_[pos++];
			subdualobj_[s] = recvbuf_[pos++];
			CoinCopyN(recvbuf_ + pos, model_->getNumSubproblemCouplingCols(s), subsolution_[s]);
			pos += model_->getNumSubproblemCouplingCols(s);
		}

		/** update problem */
		master_->updateProblem(primalBound_, dualBound_, subprimobj_, subdualobj_, subsolution_);

		/** solve problem */
		master_->solve();

		/** returns continue or stop */
		master_signal = master_->getStatus();

		/** send message */
		pos = 0;
		sendbuf_[pos] = static_cast<double>(master_signal); /**< signal */
		if (master_signal == STO_STAT_MW_CONTINUE)
		{
			pos++;
			const double * theta = master_->getSolution();
			const double * lambda = master_->getSolution() + model_->getNumSubproblems();
			for (int s = 0; s < nsubprobs; ++s)
			{
				int nlambda = model_->getNumSubproblemCouplingRows(s);
				const int * relevantRows = model_->getSubproblemCouplingRowIndices(s);
				/** theta */
				sendbuf_[pos++] = theta[subindex_[s]];
				/** lambda */
				for (int k = 0; k < nlambda; ++k)
					sendbuf_[pos++] = lambda[relevantRows[k]];
			}
		}
		MPI_Send(sendbuf_, pos, MPI_DOUBLE, status.MPI_SOURCE, MPI_ANY_TAG, comm_);
	}

	nworkers--;
	while (nworkers > 0)
	{
		/** receive message */
		MPI_Recv(recvbuf_, rcount_, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, comm_, &status);
		/** send message */
		MPI_Send(sendbuf_, 1, MPI_DOUBLE, status.MPI_SOURCE, MPI_ANY_TAG, comm_);
		/** decrement number of workers */
		nworkers--;
	}

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
}

/** run worker processes */
DSP_RTN_CODE DecDdAsyncMpi::runWorkers()
{
	if (comm_rank_ == 0)
		return STO_RTN_OK;

	BGN_TRY_CATCH

	int doContinue = 1; /**< indicate if continue to solve */

	/** loop until when the master signals stop */
	while(doContinue)
	{
		/** update subproblems */
		updateSubproblems(doContinue);

		/** Solve subproblems assigned to each process  */
		solveSubproblems(doContinue);

		/** request signal and new multiplier */
		commWithMaster(doContinue);
	}

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
}

/** update subproblems */
DSP_RTN_CODE DecDdAsyncMpi::updateSubproblems(int & doContinue)
{
	BGN_TRY_CATCH

	/** update wall clock */
	if (ticToc() < 0)
	{
		doContinue = 0;
		return STO_RTN_OK;
	}

	for (unsigned s = 0; s < subprobs_.size(); ++s)
	{
		/** set time limit */
		subprobs_[s]->setTimeLimit(CoinMin(time_remains_, parScipTimeLim_));

		/** update problems */
		subprobs_[s]->updateProblem(lambda_);
	}

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
}

/** solve subproblems */
DSP_RTN_CODE DecDdAsyncMpi::solveSubproblems(int & doContinue)
{
	BGN_TRY_CATCH

	/** update wall clock */
	if (ticToc() < 0)
	{
		doContinue = 0;
		return STO_RTN_OK;
	}

	/** Solve subproblems assigned to each process  */
	for (unsigned s = 0; s < subprobs_.size(); ++s)
	{
		bool resolve = true;
		while(resolve)
		{
			/** solve subproblem */
			subprobs_[s]->solve();

			/** solution status */
			if (ticToc() < 0 || checkStatus(subprobs_[s]) == false)
			{
				doContinue = 0;
				break;
			}

			/** set solution gap tolerance */
			if (subprobs_[s]->getPrimalBound() >= subprobs_[s]->theta_
					&& subprobs_[s]->gapTol_ > 0)
			{
				double gapTol = subprobs_[s]->gapTol_;
				gapTol = gapTol < 0.0001 ? 0.0 : gapTol * 0.1;
				subprobs_[s]->setGapTop(gapTol);
			}
			else
			{
				resolve = false;
			}
		}
	}

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
}

/** communicate with master process */
DSP_RTN_CODE DecDdAsyncMpi::commWithMaster(int & doContinue)
{
	int pos = 0;
	MPI_Status status;

	BGN_TRY_CATCH

	/** update wall clock */
	if (ticToc() < 0)
	{
		doContinue = 0;
		return STO_RTN_OK;
	}

	/** sending message:
	 * number of subproblems
	 * [for each scenario]
	 *   scenario index
	 *   primal objective
	 *   dual objective
	 *   subproblem solution for coupling columns
	 * */
	sendbuf_[pos++] = subprobs_.size();
	for (unsigned s = 0; s < subprobs_.size(); ++s)
	{
		subprobs_[s]->MPImsgbuf(sendbuf_ + pos);
		pos += 3 + subprobs_[s]->ncols_coupling_;
	}
	MPI_Send(sendbuf_, scount_, MPI_DOUBLE, 0, MPI_ANY_TAG, comm_);

	/** receive doContinue, theta and lambda */
	pos = 0;
	MPI_Recv(recvbuf_, rcount_, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, comm_, &status);
	if (static_cast<int>(recvbuf_[pos]) == STO_STAT_MW_STOP)
		doContinue = 0;
	else
	{
		pos++;
		/** continue */
		for (unsigned s = 0; s < subprobs_.size(); ++s)
		{
			subprobs_[s]->theta_ = recvbuf_[pos++];
			subprobs_[s]->updateProblem(recvbuf_ + pos);
			pos += subprobs_[s]->nrows_coupling_;
		}
	}

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
}
