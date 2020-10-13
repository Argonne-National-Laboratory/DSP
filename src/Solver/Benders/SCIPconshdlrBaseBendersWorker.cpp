/*
 * SCIPconshdlrBaseBendersWorker.cpp
 *
 *  Created on: Feb 18, 2016
 *      Author: kibaekkim
 */

// #define DSP_DEBUG
#include "Utility/DspMessage.h"
#include "Solver/Benders/SCIPconshdlrBaseBendersWorker.h"
#include "Solver/Benders/BdMW.h"

SCIPconshdlrBaseBendersWorker::SCIPconshdlrBaseBendersWorker(MPI_Comm comm)
	: comm_(comm),
	  nsubprobs_(-1),
	  recvcounts_(NULL),
	  displs_(NULL),
	  cut_index_(NULL),
	  cut_status_(NULL),
	  cutval_(NULL),
	  cutrhs_(NULL)
{
	MPI_Comm_size(comm_, &comm_size_);
	MPI_Comm_rank(comm_, &comm_rank_);
}

SCIPconshdlrBaseBendersWorker::~SCIPconshdlrBaseBendersWorker()
{
	FREE_ARRAY_PTR(recvcounts_);
	FREE_ARRAY_PTR(displs_);
	FREE_ARRAY_PTR(cut_index_);
	FREE_ARRAY_PTR(cut_status_);
	FREE_2D_ARRAY_PTR(nsubprobs_, cutval_);
	FREE_ARRAY_PTR(cutrhs_);
}

void SCIPconshdlrBaseBendersWorker::initialize(int nsubprobs)
{
	nsubprobs_ = nsubprobs;
	recvcounts_ = new int[comm_size_];
	displs_ = new int[comm_size_];
	cut_index_ = new int[nsubprobs_];
	cut_status_ = new int[nsubprobs_];
	cutval_ = new double *[nsubprobs_];
	for (int s = 0; s < nsubprobs_; ++s)
		cutval_[s] = NULL;
	cutrhs_ = new double[nsubprobs_];

	/** The master does not solve subproblem. */
	recvcounts_[0] = 0;
	displs_[0] = 0;

	/** Subproblems are assigned to each process in round-and-robin fashion. */
	int j = 0;
	for (int i = 1; i < comm_size_; ++i)
	{
		recvcounts_[i] = 0;
		for (int s = i - 1; s < nsubprobs_; s += comm_size_ - 1)
		{
			recvcounts_[i]++;
			cut_index_[j++] = s;
		}
		displs_[i] = i == 0 ? 0 : displs_[i - 1] + recvcounts_[i - 1];
		DSPdebugMessage("recvcounts_[%d] %d displs_[%d] %d\n", i, recvcounts_[i], i, displs_[i]);
	}
	DSPdebug(DspMessage::printArray(nsubprobs_, cut_index_));
}

void SCIPconshdlrBaseBendersWorker::generateCutsBase(
	int nsubprobs,			   /**< [in] number of subproblems */
	int nvars,				   /**< [in] number of variables */
	int naux,				   /**< [in] number of auxiliary variables */
	double *x,				   /**< [in] master solution */
	const double *probability, /**< [in] probability */
	OsiCuts *cuts /**< [out] cuts generated */)
{
	OsiCuts cuts_collected;

	BGN_TRY_CATCH

	/** Tell workers to generate cuts */
	int message = BdMW::MASTER_NEEDS_CUTS;
	MPI_Bcast(&message, 1, MPI_INT, 0, comm_);

	/** Send solutions to the workers */
	MPI_Bcast(x, nvars, MPI_DOUBLE, 0, comm_);
	// printf("x:\n");
	// DspMessage::printArray(nvars, x);

	/** Collect cut generation stutus */
	MPI_Gatherv(NULL, 0, MPI_INT, cut_status_, recvcounts_, displs_, MPI_INT, 0, comm_);

	/** Collect cuts */
	MPIgatherOsiCuts(comm_, *cuts, cuts_collected);
	DSPdebugMessage("[%d]: Collected %d cuts\n", comm_rank_, cuts_collected.sizeCuts());

	for (int i = 0; i < cuts_collected.sizeCuts(); ++i)
	{
		FREE_ARRAY_PTR(cutval_[i]);
		OsiRowCut *rc = cuts_collected.rowCutPtr(i);
		cutval_[i] = rc->row().denseVector(nvars);
		cutrhs_[i] = rc->lb();
	}

	/** aggregate cuts */
	aggregateCutsBase(nsubprobs, nvars, naux, probability, cuts);
	DSPdebug(cuts->printCuts());
	// cuts->printCuts();

	END_TRY_CATCH(;)

#undef FREE_MEMORY
}

void SCIPconshdlrBaseBendersWorker::aggregateCutsBase(
	int nsubprobs,			   /**< [in] number of subproblems */
	int nvars,				   /**< [in] number of variables */
	int naux,				   /**< [in] number of auxiliary variables */
	const double *probability, /**< [in] probability */
	OsiCuts *cuts /**< [out] cuts generated */)
{
#define FREE_MEMORY                 \
	FREE_2D_ARRAY_PTR(naux, aggval) \
	FREE_ARRAY_PTR(aggrhs)

	bool isInfeasible = false; /**< indicating whether there is primal infeasibility or not */
	double **aggval = NULL;	   /** aggregated dense cut coefficients */
	double *aggrhs = NULL;	   /** aggregated cut rhs */
	CoinPackedVector vec;

	BGN_TRY_CATCH

	/** allocate memory */
	aggval = new double *[naux];
	for (int s = naux - 1; s >= 0; --s)
	{
		aggval[s] = new double[nvars];
		CoinZeroN(aggval[s], nvars);
	}
	aggrhs = new double[naux];
	CoinZeroN(aggrhs, naux);

	for (int i = nsubprobs - 1; i >= 0; --i)
	{
		/** generate feasibility cut */
		if (cut_status_[i] == DSP_STAT_PRIM_INFEASIBLE)
		{
			/** set cut body */
			for (int j = 0; j < nvars; ++j)
				if (fabs(cutval_[i][j]) > 1.0e-8)
					vec.insert(j, cutval_[i][j]);

			OsiRowCut fcut;
			fcut.setRow(vec);
			fcut.setUb(COIN_DBL_MAX);
			fcut.setLb(cutrhs_[i]);

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

		int ind_aux = cut_index_[i] % naux;

		/** calculate weighted aggregation of cuts */
		for (int j = 0; j < nvars; ++j)
			aggval[ind_aux][j] += cutval_[i][j] * probability[cut_index_[i]];
		aggrhs[ind_aux] += cutrhs_[i] * probability[cut_index_[i]];
	}

	/** We generate optimality cuts only if there is no feasibility cut generated. */
	if (isInfeasible == false)
	{
		/** construct cuts to pass */
		for (int s = 0; s < naux; ++s)
		{
			/** auxiliary variable coefficient */
			aggval[s][nvars - naux + s] = 1;

			/** initialize vector */
			vec.clear();

			/** set it as sparse */
			for (int j = 0; j < nvars; ++j)
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
