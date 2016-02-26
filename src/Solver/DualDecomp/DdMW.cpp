/*
 * DdMW.cpp
 *
 *  Created on: Feb 9, 2016
 *      Author: kibaekkim
 */

#include "Solver/DualDecomp/DdMW.h"

DdMW::DdMW(
		MPI_Comm comm,
		Master * master,
		Worker * worker):
	BaseMasterWorker(comm),
	master_(master), worker_(worker) {}

DdMW::~DdMW() {}

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
	int scount = 0;

	BGN_TRY_CATCH

	/** buffer size */
	scount_ = comm_rank_ == 0 ? master_->getSendCount() : worker_->getSendCount();
	rcount_ = comm_rank_ == 0 ? master_->getRecvCount() : worker_->getRecvCount();

	if (sync_)
	{
		if (comm_rank_ == 0)
		{
			rcounts_ = new int [comm_size_];
			sdispls_ = new int [comm_size_];
			rdispls_ = new int [comm_size_];
			scount_ = 0;
			rcount_ = 0;
		}
		MPI_Gather(&scount_, 1, MPI_INT, rcounts_, comm_size_, MPI_INT, 0, comm_);
		MPI_Gather(&rcount_, 1, MPI_INT, scounts_, comm_size_, MPI_INT, 0, comm_);
		sdispls_[0] = 0;
		rdispls_[0] = 0;
		for (int i = 1; i < comm_size_; ++i)
		{
			sdispls_[i] = sdispls_[i-1] + scounts_[i-1];
			rdispls_[i] = rdispls_[i-1] + rcounts_[i-1];
		}
	}
	else
	{
		/** allocate memory */
		recvbuf_ = new double [rcount_];
	}

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
}

/** run master process */
STO_RTN_CODE DdMW::runMaster()
{
	if (comm_rank_ != 0)
		return STO_RTN_OK;

	MPI_Status status;
	STO_RTN_CODE signal = STO_STAT_MW_CONTINUE;
	int nworkers = comm_size_ - 1; /**< number of workers */
	double garbage = -1;

	BGN_TRY_CATCH

	while (signal == STO_STAT_MW_CONTINUE)
	{
		/** receive message */
		if (sync_)
		{
			STO_RTN_CODE ssignal = signal;
			MPI_Reduce(&ssignal, &signal, 1, MPI_INT, MPI_MIN, 0, comm_);
			MPI_Gatherv(NULL, 0, MPI_DOUBLE, recvbuf_, rcounts_, rdispls_, MPI_DOUBLE, 0, comm_);
		}
		else
		{
			MPI_Recv(recvbuf_, rcount_, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, comm_, &status);
			signal = status.MPI_TAG;
		}
		master_->recvMessage(status.MPI_SOURCE, status.count, recvbuf_);

		if (signal == STO_STAT_MW_STOP)
			break;

		/** update problem */
		master_->updateProblem();

		/** solve problem */
		master_->solve();

		/** returns continue or stop signal */
		signal = master_->getStatus();

		/** send message */
		if (sync_)
		{
			MPI_Bcast(&signal, 1, MPI_INT, 0, comm_);
			if (signal == STO_STAT_MW_CONTINUE)
				MPI_Scatterv(master_->getSendMessage(), scounts_, sdispls_,
						MPI_DOUBLE, NULL, 0, MPI_DOUBLE, 0, comm_);
		}
		else
			MPI_Send(master_->getSendMessage(), master_->getSendCount(),
					MPI_DOUBLE, status.MPI_SOURCE, signal, comm_);
	}

	/** sending stop signals one-by-one for asynchronous */
	if (!sync_)
	{
		nworkers--;

		/** sending stop signals to all the worker processes */
		while (nworkers > 0)
		{
			/** receive message */
			MPI_Recv(recvbuf_, rcount_, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, comm_, &status);
			/** send garbage with stop signal as a tag */
			MPI_Send(&garbage, 1, MPI_DOUBLE, status.MPI_SOURCE, signal, comm_);
			/** decrement number of workers */
			nworkers--;
		}
	}

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
}

/** ruw worker processes */
STO_RTN_CODE DdMW::runWorker()
{
	if (comm_rank_ == 0)
		return STO_RTN_OK;

	MPI_Status status;
	STO_RTN_CODE signal = STO_STAT_MW_CONTINUE;
	int doContinue = 1; /**< indicate if continue to solve */

	BGN_TRY_CATCH

	/** loop until when the master signals stop */
	while(doContinue)
	{
		/** Solve subproblems assigned to each process  */
		worker_->solve();

		if (sync_)
		{
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

			/** send message to the master */
			int rsignal;
			MPI_Reduce(&signal, &rsignal, 1, MPI_INT, MPI_MIN, 0, comm_);
			MPI_Gatherv(worker_->getSendMessage(), worker_->getSendCount(),
					MPI_DOUBLE, NULL, NULL, NULL, MPI_DOUBLE, 0, comm_);

			/** receive message from the master */
			MPI_Bcast(&signal, 1, MPI_INT, 0, comm_);
			if (signal == STO_STAT_MW_CONTINUE)
			{
				MPI_Scatterv(NULL, NULL, NULL, MPI_DOUBLE,
						recvbuf_, rcount_, MPI_DOUBLE, 0, comm_);
			}
		}
		else
		{
			/** worker status */
			if (worker_->getStatus() == STO_STAT_MW_STOP)
			{
				/** send stop signal to master */
				MPI_Send(NULL, 0, MPI_DOUBLE, 0, STO_STAT_MW_STOP, comm_);
				/** stop */
				doContinue = 0;
				break;
			}

			/** send message to the master */
			MPI_Send(worker_->getSendMessage(), worker_->getSendCount(),
					MPI_DOUBLE, 0, MPI_ANY_TAG, comm_);

			/** receive message from the master */
			MPI_Recv(recvbuf_, rcount_, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, comm_, &status);
			signal = status.MPI_TAG;
		}
		if (signal == STO_STAT_MW_STOP)
			doContinue = 0;
		else
		{
			/** update subproblems */
			worker_->recvMessage(status.count, recvbuf_);
		}
	}

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
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

