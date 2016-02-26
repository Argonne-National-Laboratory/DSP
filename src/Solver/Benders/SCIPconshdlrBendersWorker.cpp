/*
 * SCIPconshdlrBendersWorker.cpp
 *
 *  Created on: Feb 18, 2016
 *      Author: kibaekkim
 */

#include "Utility/StoMessage.h"
#include "Solver/Benders/SCIPconshdlrBendersWorker.h"
#include "Solver/Benders/BdMW.h"

SCIPconshdlrBendersWorker::SCIPconshdlrBendersWorker(SCIP * scip, int sepapriority, MPI_Comm comm):
	SCIPconshdlrBenders(scip, sepapriority),
	comm_(comm), recvcounts_(NULL), displs_(NULL),
	cut_indices_(NULL), cut_status_(NULL)
{
	MPI_Comm_size(comm_, &comm_size_);
	MPI_Comm_rank(comm_, &comm_rank_);
}

SCIPconshdlrBendersWorker::~SCIPconshdlrBendersWorker() {}

/** destructor of constraint handler to free user data (called when SCIP is exiting) */
SCIP_DECL_CONSFREE(SCIPconshdlrBendersWorker::scip_free)
{
	for (int j = 0; j < nvars_; ++j)
		vars_[j] = NULL;
	SCIPfreeMemoryArray(scip_, &vars_);
	nvars_ = 0;
	naux_ = 0;

	FREE_ARRAY_PTR(recvcounts_);
	FREE_ARRAY_PTR(displs_);
	FREE_ARRAY_PTR(cut_indices_);
	FREE_ARRAY_PTR(cut_status_);

	return SCIP_OKAY;
}

/** clone method which will be used to copy constraint handler and variable pricer objects */
SCIP_DECL_CONSHDLRCLONE(scip::ObjProbCloneable* SCIPconshdlrBendersWorker::clone)
{
	*valid = true;
	SCIPconshdlrBendersWorker * conshdlrclone = new SCIPconshdlrBendersWorker(scip, scip_sepapriority_, comm_);
	conshdlrclone->setOriginalVariables(nvars_, vars_, naux_);
	return conshdlrclone;
}

void SCIPconshdlrBendersWorker::setNumSubprobs(int nsubprobs)
{
	nsubprobs_   = nsubprobs;
	recvcounts_  = new int [comm_size_];
	displs_      = new int [comm_size_];
	cut_indices_ = new int [nsubprobs_];
	cut_status_  = new int [nsubprobs_];

	/** The master does not solve subproblem. */
	recvcounts_[0] = 0;
	displs_[0] = 0;

	/** Subproblems are assigned to each process in round-and-robin fashion. */
	for (int i = 1, j = 0; i < comm_size_; ++i)
	{
		recvcounts_[i] = 0;
		for (int s = i; s < nsubprobs_; s += comm_size_-1)
		{
			recvcounts_[i]++;
			cut_indices_[j++] = s;
		}
		displs_[i] = i == 0 ? 0 : displs_[i-1] + recvcounts_[i-1];
		DSPdebugMessage("recvcounts_[%d] %d displs_[%d] %d\n", i, recvcounts_[i], i, displs_[i]);
	}
}

void SCIPconshdlrBendersWorker::generateCuts(
		int size,      /**< [in] size of x */
		double * x,    /**< [in] master solution */
		int where,     /**< [in] where to be called */
		OsiCuts * cuts /**< [out] cuts generated */)
{
#define FREE_MEMORY                      \
	FREE_2D_ARRAY_PTR(nsubprobs, cutval) \
	FREE_ARRAY_PTR(cutrhs)

	OsiCuts cuts_collected;
	int nsubprobs = bdsub_->getNumSubprobs();
	double ** cutval = NULL;   /** dense cut coefficients for each subproblem */
	double *  cutrhs = NULL;   /** cut rhs for each subproblem */

	BGN_TRY_CATCH

	/** Tell workers to generate cuts */
	int message = BdMW::MASTER_NEEDS_CUTS;
	MPI_Bcast(&message, 1, MPI_INT, 0, comm_);

	/** Send solutions to the workers */
	MPI_Bcast(x, nvars_, MPI_DOUBLE, 0, comm_);

	/** Collect cut generation stutus */
	MPI_Gatherv(NULL, 0, MPI_INT, cut_status_, recvcounts_, displs_, MPI_INT, 0, comm_);

	/** Collect cuts */
	MPIgatherOsiCuts(comm_, *cuts, cuts_collected);
	DSPdebugMessage("[%d]: Collected %d cuts\n", comm_rank_, cuts_collected.sizeCuts());

	/** allocate memory */
	cutval = new double * [nsubprobs];
	cutrhs = new double [nsubprobs];

	for (int i = 0; i < cuts_collected.sizeCuts(); ++i)
	{
		OsiRowCut * rc = cuts_collected.rowCutPtr(i);
		cutval[i] = rc->row().denseVector(nvars_);
		cutrhs[i] = rc->lb();
	}

	/** aggregate cuts */
	aggregateCuts(cutval, cutrhs, cuts);
	DSPdebug(cuts->printCuts());

	END_TRY_CATCH(FREE_MEMORY)

	FREE_MEMORY

#undef FREE_MEMORY
}
