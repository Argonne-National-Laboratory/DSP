// #define DSP_DEBUG
#include "Model/TssModel.h"
#include "Solver/DualDecomp/DdDroWorkerUBMpi.h"

DdDroWorkerUBMpi::DdDroWorkerUBMpi(
	MPI_Comm comm,	 /**< MPI communicator */
	DecModel *model, /**< model pointer */
	DspParams *par,	 /**< parameter pointer */
	DspMessage *message /**< message pointer */)
	: DdDroWorkerUB(model, par, message), comm_(comm)
{
	MPI_Comm_size(comm_, &comm_size_);
	MPI_Comm_rank(comm_, &comm_rank_);
	DSPdebugMessage("comm_size_ %d comm_rank_ %d\n", comm_size_, comm_rank_);
}

DSP_RTN_CODE DdDroWorkerUBMpi::init()
{
	BGN_TRY_CATCH
	status_ = DSP_STAT_MW_CONTINUE;

	/** Create the stochastic upper bounding subproblems.
	 * We still need to solve these subproblems,
	 * in addition to the DRO upper bounding problem.
	 */
	DSP_RTN_CHECK_THROW(DdWorkerUB::createProblem());

	/** Create the DRO upper bounding problem. 
	 * This problem is solved by the root process only.
	*/
	if (comm_rank_ == 0)
		DSP_RTN_CHECK_THROW(createProblem());

	DSP_RTN_CHECK_THROW(setObjective());

	END_TRY_CATCH_RTN(;, DSP_RTN_ERR)
	return DSP_RTN_OK;
}

DSP_RTN_CODE DdDroWorkerUBMpi::solve()
{
#define FREE_MEMORY            \
	FREE_ARRAY_PTR(sendbuf)    \
	FREE_ARRAY_PTR(recvbuf)    \
	FREE_ARRAY_PTR(recvcounts) \
	FREE_ARRAY_PTR(displs)     \
	FREE_ARRAY_PTR(subprob_ids)

	double *sendbuf = NULL;
	double *recvbuf = NULL;
	int *recvcounts = NULL;
	int *displs = NULL;
	int *subprob_ids = NULL;

	double cputime;
	double walltime;

	BGN_TRY_CATCH

	double primobj = 0.0;
	double dualobj = 0.0;
	double total_cputime = 0.0;
	double total_walltime = 0.0;
	int nsubprobs = par_->getIntPtrParamSize("ARR_PROC_IDX");

	for (unsigned s = 0; s < nsubprobs; ++s)
	{
		cputime = CoinCpuTime();
		walltime = CoinGetTimeOfDay();

		/** set time limit */
		osi_[s]->setTimeLimit(
			CoinMin(CoinMax(0.01, time_remains_),
					par_->getDblParam("DD/SUB/TIME_LIM")));

		/** solve */
		osi_[s]->solve();

		/** check status. there might be unexpected results. */
		int status = osi_[s]->status();
		DSPdebugMessage("status = %d\n", status);
		switch (status)
		{
		case DSP_STAT_OPTIMAL:
		case DSP_STAT_LIM_ITERorTIME:
		case DSP_STAT_STOPPED_GAP:
		case DSP_STAT_STOPPED_NODE:
		case DSP_STAT_STOPPED_TIME:
			break;
		default:
			status_ = DSP_STAT_MW_STOP;
			message_->print(10,
							"Warning: subproblem %d solution status is %d\n", s,
							status);
			break;
		}
		if (status_ == DSP_STAT_MW_STOP)
		{
			primobj = COIN_DBL_MAX;
			dualobj = -COIN_DBL_MAX;
			break;
		}

		CoinCopyN(osi_[s]->si_->getColSolution(), osi_[s]->si_->getNumCols(), &primsols_[s][0]);
		total_cputime += CoinCpuTime() - cputime;
		total_walltime += CoinGetTimeOfDay() - walltime;

		/** consume time */
		time_remains_ -= CoinGetTimeOfDay() - walltime;
	}

	cputime = CoinCpuTime();
	walltime = CoinGetTimeOfDay();

	/** solve the DRO UB problem */
	if (status_ != DSP_STAT_MW_STOP)
	{
		/** memory allocation for MPI communication */
		sendbuf = new double[nsubprobs];
		if (comm_rank_ == 0)
		{
			recvbuf = new double[model_->getNumSubproblems()];
			recvcounts = new int[comm_size_];
			displs = new int[comm_size_];
			subprob_ids = new int[model_->getNumSubproblems()];
		}

		/** recvcounts stores the number of subproblems for each process. */
		// DSPdebugMessage("Rank %d: nsubprobs %d\n", comm_rank_, nsubprobs);
		MPI_Gather(&nsubprobs, 1, MPI_INT, recvcounts, 1, MPI_INT, 0, comm_);
		if (comm_rank_ == 0)
		{
			for (int i = 0; i < comm_size_; ++i)
			{
				displs[i] = i == 0 ? 0 : displs[i - 1] + recvcounts[i - 1];
			}
			DSPdebugMessage("recvcounts:\n");
			DSPdebug(message_->printArray(comm_size_, recvcounts));
			DSPdebugMessage("displs:\n");
			DSPdebug(message_->printArray(comm_size_, displs));
		}

		/** assign send buffer */
		for (int i = 0; i < nsubprobs; ++i)
			sendbuf[i] = osi_[i]->si_->getObjValue();
		// DSPdebugMessage("Rank %d: sendbuf:\n", comm_rank_);
		// DSPdebug(message_->printArray(nsubprobs, sendbuf));

		/** Gather subproblem indices */
		MPI_Gatherv(par_->getIntPtrParam("ARR_PROC_IDX"), nsubprobs, MPI_INT,
					subprob_ids, recvcounts, displs, MPI_INT, 0, comm_);

		/** Gather subproblem objective values */
		MPI_Gatherv(sendbuf, nsubprobs, MPI_DOUBLE, recvbuf, recvcounts, displs, MPI_DOUBLE, 0, comm_);

		if (comm_rank_ == 0)
		{
			// DSPdebugMessage("subproblem IDs:\n");
			// DSPdebug(message_->printArray(model_->getNumSubproblems(), subprob_ids));
			// DSPdebugMessage("subproblem objs:\n");
			// DSPdebug(message_->printArray(model_->getNumSubproblems(), recvbuf));

			/** Set right-hand sides for the upper bounding problem */
			for (int k = 0; k < model_->getNumSubproblems(); ++k)
			{
				DSPdebugMessage("subprob_ids[%d] = %d, rlbd[%d] = %e\n", k, subprob_ids[k],
								subprob_ids[k] * model_->getNumReferences(),
								recvbuf[k]);
				for (int s = 0; s < model_->getNumReferences(); ++s)
				{
					osi_dro_->si_->setRowLower(
						subprob_ids[k] * model_->getNumReferences() + s,
						recvbuf[k]);
				}
			}
			osi_dro_->solve();

			if (osi_dro_->si_->isProvenOptimal())
			{
				primobj = osi_dro_->getPrimObjValue();
				dualobj = osi_dro_->getPrimObjValue();
			}
			else
			{
				primobj = COIN_DBL_MAX;
				dualobj = -COIN_DBL_MAX;
			}
		}
		else

		{
			primobj = 0.0;
			dualobj = 0.0;
		}
	}

	/** get primal objective */
	ub_ = primobj;
	DSPdebugMessage("ub_ = %e\n", ub_);
	DSPdebugMessage("status_ %d\n", status_);

	MPI_Barrier(comm_);

	total_cputime += CoinCpuTime() - cputime;
	total_walltime += CoinGetTimeOfDay() - walltime;
	time_remains_ -= CoinGetTimeOfDay() - walltime;

	/** update statistics */
	s_statuses_.push_back(status_);
	s_primobjs_.push_back(primobj);
	s_dualobjs_.push_back(dualobj);
	s_cputimes_.push_back(total_cputime);
	s_walltimes_.push_back(total_walltime);

	END_TRY_CATCH_RTN(FREE_MEMORY, DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}
