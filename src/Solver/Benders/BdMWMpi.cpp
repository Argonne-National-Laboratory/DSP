/*
 * BdMWMpi.cpp
 *
 *  Created on: Jul 10, 2016
 *      Author: kibaekkim
 */

#define DSP_DEBUG
#include "Solver/Benders/BdMWMpi.h"
#include "Solver/Benders/SCIPconshdlrBendersWorker.h"
#include "SolverInterface/SolverInterfaceScip.h"

BdMWMpi::BdMWMpi(
		MPI_Comm comm,
		DecModel* model, /**< model pointer */
		DspParams* par, /**< parameters */
		DspMessage* message):
BdMW(model, par, message),
comm_(comm)
{
	MPI_Comm_rank(comm, &comm_rank_);
	MPI_Comm_size(comm, &comm_size_);
}

BdMWMpi::~BdMWMpi() {
	// TODO Auto-generated destructor stub
}

DSP_RTN_CODE BdMWMpi::init()
{
	BGN_TRY_CATCH

	if (comm_rank_ == 0)
	{
		/** create and initialize master */
		master_ = new BdMaster(par_, model_, message_);
		DSP_RTN_CHECK_THROW(master_->init());
	}
	else
	{
		/** set parameters */
		vector<int> subprob_indices;
		distIndices(model_->getNumSubproblems(), comm_size_ - 1, comm_rank_, 1, subprob_indices);
		par_->setIntPtrParamSize("ARR_PROC_IDX", subprob_indices.size());
		for (unsigned s = 0; s < subprob_indices.size(); ++s)
			par_->setIntPtrParam("ARR_PROC_IDX", s, subprob_indices[s]);

		/** create and initialize worker */
		worker_ = new BdWorker(par_, model_, message_);
		DSP_RTN_CHECK_THROW(worker_->init());
	}
#ifdef DSP_DEBUG
	for (int i = 1; i < comm_size_; ++i)
	{
		if (i == comm_rank_)
		{
			DSPdebugMessage("Rank %d: ARR_PROC_IDX:\n", comm_rank_);
			DspMessage::printArray(par_->getIntPtrParamSize("ARR_PROC_IDX"),
					par_->getIntPtrParam("ARR_PROC_IDX"));
		}
		MPI_Barrier(comm_);
	}
#endif

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE BdMWMpi::finalize()
{
	BGN_TRY_CATCH

	if (comm_rank_ == 0)
	{
		master_->finalize();
		FREE_PTR(master_);
	}
	else
	{
		worker_->finalize();
		FREE_PTR(worker_);
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE BdMWMpi::run()
{
	BGN_TRY_CATCH

	/** run master process */
	DSP_RTN_CHECK_THROW(runMaster());

	/** run worker processes */
	DSP_RTN_CHECK_THROW(runWorker());

	DSPdebugMessage("Rank %d: End of run()\n", comm_rank_);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE BdMWMpi::runMaster()
{
	if (!master_)
		return DSP_RTN_OK;

	BGN_TRY_CATCH

	DSPdebugMessage("Rank %d: Run Master\n", comm_rank_);
	/** set constraint handler */
	master_->setConshdlr(constraintHandler());

	/** solve */
	DSP_RTN_CHECK_THROW(master_->solve());

	/** Tell workers we are done */
	int message = MASTER_STOPPED;
	MPI_Bcast(&message, 1, MPI_INT, 0, comm_);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

SCIPconshdlrBenders* BdMWMpi::constraintHandler()
{
	SCIPconshdlrBenders * conshdlr = NULL;

	BGN_TRY_CATCH

	/** get solver interface */
	SolverInterfaceScip * si = dynamic_cast<SolverInterfaceScip*>(master_->getSiPtr());

	/** MPI Benders */
	conshdlr = new SCIPconshdlrBendersWorker(si->getSCIP(), par_->getIntParam("BD/CUT_PRIORITY"), comm_);
	conshdlr->setDecModel(model_);
	conshdlr->setOriginalVariables(si->getNumCols(), si->getSCIPvars(), par_->getIntParam("BD/NUM_CUTS_PER_ITER"));

	END_TRY_CATCH_RTN(;,NULL)

	return conshdlr;
}

DSP_RTN_CODE BdMWMpi::runWorker()
{
#define FREE_MEMORY \
	FREE_ARRAY_PTR(solution); \
	FREE_2D_ARRAY_PTR(parProcIdxSize,cutval); \
	FREE_ARRAY_PTR(status); \
	FREE_ARRAY_PTR(cutrhs);

	if (!worker_)
		return DSP_RTN_OK;

	DSPdebugMessage("Beginning of runWorker()\n");

	OsiCuts cuts, tempcuts;
	int message;
	int ncols;
	double * solution = NULL;
	double ** cutval = NULL;
	double * cutrhs = NULL;
	CoinPackedVector vec;
	int * status = NULL;

	/** internal pointer */
	BdSub * bdsub = worker_->getBdSubPtr();
	int parProcIdxSize = par_->getIntPtrParamSize("ARR_PROC_IDX");

	BGN_TRY_CATCH

	DSPdebugMessage("Rank %d: Run Worker\n", comm_rank_);

	ncols = model_->getNumSubproblemCouplingCols(0) + par_->getIntParam("BD/NUM_CUTS_PER_ITER");
	solution = new double [ncols];
	cutval = new double * [parProcIdxSize];
	cutrhs = new double [parProcIdxSize];
	status = new int [parProcIdxSize];
	for (int i = 0; i < parProcIdxSize; ++i)
		cutval[i] = NULL;

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
		for (int s = 0; s < bdsub->getNumSubprobs(); ++s)
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
			status[s] = bdsub->getStatus(s);
			DSPdebugMessage("[%d]: status[%d] %d\n", comm_rank_, s, status[s]);
		}
		DSPdebugMessage("[%d]: Found %d cuts\n", comm_rank_, cuts.sizeCuts());

		/** Send cut generation status to the master */
		DSPdebugMessage("parProcIdxSize %d\n", parProcIdxSize);
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

	DSPdebugMessage("Rank %d: End of runWorker()\n", comm_rank_);

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	DSPdebugMessage("Rank %d: Free memory\n", comm_rank_);

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}
