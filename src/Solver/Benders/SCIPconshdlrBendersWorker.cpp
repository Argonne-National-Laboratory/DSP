/*
 * SCIPconshdlrBendersWorker.cpp
 *
 *  Created on: Feb 18, 2016
 *      Author: kibaekkim
 */

// #define DSP_DEBUG
#include "Utility/DspMessage.h"
#include "Solver/Benders/SCIPconshdlrBendersWorker.h"
#include "Solver/Benders/BdMW.h"

SCIPconshdlrBendersWorker::SCIPconshdlrBendersWorker(
		SCIP * scip,
		int sepapriority,
		bool is_integral_recourse,
		MPI_Comm comm):
SCIPconshdlrBenders(scip, "Benders", sepapriority),
comm_(comm),
recvcounts_(NULL),
displs_(NULL),
cut_index_(NULL),
cut_status_(NULL),
is_integral_recourse_(is_integral_recourse)
{
	DSPdebugMessage("is_integral_recourse_? %s\n", is_integral_recourse_ ? "yes" : "no");
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
	FREE_ARRAY_PTR(cut_index_);
	FREE_ARRAY_PTR(cut_status_);

	return SCIP_OKAY;
}

void SCIPconshdlrBendersWorker::setDecModel(DecModel * model)
{
	model_        = model;
	int nsubprobs = model_->getNumSubproblems();
	recvcounts_   = new int [comm_size_];
	displs_       = new int [comm_size_];
	cut_index_    = new int [nsubprobs];
	cut_status_   = new int [nsubprobs];

	/** The master does not solve subproblem. */
	recvcounts_[0] = 0;
	displs_[0] = 0;

	/** Subproblems are assigned to each process in round-and-robin fashion. */
	for (int i = 1; i < comm_size_; ++i)
	{
		recvcounts_[i] = 0;
		for (int s = i - 1; s < nsubprobs; s += comm_size_-1)
		{
			recvcounts_[i]++;
			cut_index_[s] = i;
		}
		displs_[i] = i == 0 ? 0 : displs_[i-1] + recvcounts_[i-1];
		DSPdebugMessage("recvcounts_[%d] %d displs_[%d] %d\n", i, recvcounts_[i], i, displs_[i]);
		DSPdebug(DspMessage::printArray(nsubprobs, cut_index_));
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
	double ** cutval = NULL;   /** dense cut coefficients for each subproblem */
	double *  cutrhs = NULL;   /** cut rhs for each subproblem */
	int nsubprobs = model_->getNumSubproblems();

	BGN_TRY_CATCH

	/** Tell workers to generate cuts */
	int message = BdMW::MASTER_NEEDS_CUTS;
	MPI_Bcast(&message, 1, MPI_INT, 0, comm_);

	/** Send solutions to the workers */
	MPI_Bcast(x, nvars_, MPI_DOUBLE, 0, comm_);
	// printf("x:\n");
	// DspMessage::printArray(nvars_, x);

	/** Collect cut generation stutus */
	MPI_Gatherv(NULL, 0, MPI_INT, cut_status_, recvcounts_, displs_, MPI_INT, 0, comm_);

	/** Collect cuts */
	MPIgatherOsiCuts(comm_, *cuts, cuts_collected);
	DSPdebugMessage("[%d]: Collected %d cuts\n", comm_rank_, cuts_collected.sizeCuts());

	/** TODO: Collect integer feasibility status to the master */

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
	// cuts->printCuts();

	END_TRY_CATCH(FREE_MEMORY)

	FREE_MEMORY

#undef FREE_MEMORY
}

void SCIPconshdlrBendersWorker::aggregateCuts(
		double ** cutvec, /**< [in] cut vector */
		double *  cutrhs, /**< [in] cut right-hand side */
		OsiCuts * cuts    /**< [out] cuts generated */)
{
#define FREE_MEMORY                      \
	FREE_2D_ARRAY_PTR(naux_, aggval)     \
	FREE_ARRAY_PTR(aggrhs)

	int nsubprobs = model_->getNumSubproblems();
	bool isInfeasible = false; /**< indicating whether there is primal infeasibility or not */
	double ** aggval = NULL;   /** aggregated dense cut coefficients */
	double *  aggrhs = NULL;   /** aggregated cut rhs */
	CoinPackedVector vec;

	BGN_TRY_CATCH

	/** allocate memory */
	aggval = new double * [naux_];
	for (int s = naux_ - 1; s >= 0; --s)
	{
		aggval[s] = new double [nvars_];
		CoinZeroN(aggval[s], nvars_);
	}
	aggrhs = new double [naux_];
	CoinZeroN(aggrhs, naux_);

	for (int i = nsubprobs - 1; i >= 0; --i)
	{
		/** generate feasibility cut */
		if (cut_status_[i] == DSP_STAT_PRIM_INFEASIBLE)
		{
			/** set cut body */
			for (int j = 0; j < nvars_; ++j)
				if (fabs(cutvec[i][j]) > 1.0e-8)
					vec.insert(j, cutvec[i][j]);

			OsiRowCut fcut;
			fcut.setRow(vec);
			fcut.setUb(COIN_DBL_MAX);
			fcut.setLb(cutrhs[i]);

			cuts->insert(fcut);
			isInfeasible = true;
			break;
		}
		
		/** When some subproblems were primal infeasible, the rest are not solved. Then, just skip them. */
		if (cut_status_[i] == DSP_STAT_NOT_SOLVED)
			break;

		if (cut_status_[i] != DSP_STAT_OPTIMAL)
		{
			printf("Error: Subproblem %d returns unexpected status %d\n",
					cut_index_[i], cut_status_[i]);
			for (int j = 0; j < cuts->sizeCuts(); ++j)
				delete cuts->rowCutPtr(i);
			cuts->dumpCuts();
			break;
		}

		int ind_aux = cut_index_[i] % naux_;

		/** calculate weighted aggregation of cuts */
		for (int j = 0; j < nvars_; ++j)
			aggval[ind_aux][j] += cutvec[i][j];
		aggrhs[ind_aux] += cutrhs[i];
	}

	/** We generate optimality cuts only if there is no feasibility cut generated. */
	if (isInfeasible == false)
	{
		/** construct cuts to pass */
		for (int s = 0; s < naux_; ++s)
		{
			/** auxiliary variable coefficient */
			aggval[s][nvars_ - naux_ + s] = 1;

			/** initialize vector */
			vec.clear();

			/** set it as sparse */
			for (int j = 0; j < nvars_; ++j)
				if (fabs(aggval[s][j]) > 1e-10)
					vec.insert(j, aggval[s][j]);

			if (fabs(aggrhs[s]) < 1E-10)
				aggrhs[s] = 0.0;

			OsiRowCut rc;
			rc.setRow(vec);
			rc.setUb(COIN_DBL_MAX); /** TODO: for minimization */
			rc.setLb(aggrhs[s]);

			cuts->insert(rc);
		}
	}

	END_TRY_CATCH(FREE_MEMORY)

	FREE_MEMORY

#undef FREE_MEMORY
}

SCIP_RETCODE SCIPconshdlrBendersWorker::evaluateRecourse(
		SCIP * scip,    /**< [in] scip pointer */
		SCIP_SOL * sol, /**< [in] solution to evaluate */
		double * values /**< [out] evaluated recourse values */) {
	
	/**< current solution */
	SCIP_Real * vals = NULL;

	/** recourse values collected from ranks but not sorted */
	double* recourse = NULL;

	/** Tell workers to evaluate recourse */
	int message = BdMW::MASTER_EVALUATES_RECOURSE;
	MPI_Bcast(&message, 1, MPI_INT, 0, comm_);

	/** allocate memory */
	SCIP_CALL(SCIPallocMemoryArray(scip, &vals, nvars_));

	/** get current solution */
	SCIP_CALL(SCIPgetSolVals(scip, sol, nvars_, vars_, vals));

	/** Send solutions to the workers */
	MPI_Bcast(vals, nvars_, MPI_DOUBLE, 0, comm_);

	/** check return value */
	int ret = 0, allret = 0;
	MPI_Allreduce(&ret, &allret, 1, MPI_INT, MPI_SUM, comm_);
	if (allret > 0) return SCIP_ERROR;

	/** receive recourse evaluations */
	recourse = new double [model_->getNumSubproblems()];
	MPI_Gatherv(NULL, 0, MPI_DOUBLE, recourse, recvcounts_, displs_, MPI_DOUBLE, 0, comm_);

	/** sort by index */
	int j = 0;
	for (int i = 1; i < comm_size_; ++i) {
		for (int s = i - 1; s < model_->getNumSubproblems(); s += comm_size_-1) {
			DSPdebugMessage("values[%d] = %e\n", s, recourse[j]);
			values[s] = recourse[j++];
		}
	}

	FREE_ARRAY_PTR(recourse)

	/** free memory */
	SCIPfreeMemoryArray(scip, &vals);

	return SCIP_OKAY;
}
