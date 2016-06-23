/*
 * SCIPconshdlrBendersWorker.cpp
 *
 *  Created on: Feb 18, 2016
 *      Author: kibaekkim
 */

#include "SolverInterface/SCIPconshdlrBendersWorker.h"
#include <Utility/DspMessage.h>

/** destructor of constraint handler to free user data (called when SCIP is exiting) */
SCIP_DECL_CONSFREE(SCIPconshdlrBendersWorker::scip_free)
{
	for (int j = 0; j < nvars_; ++j)
		vars_[j] = NULL;
	SCIPfreeMemoryArray(scip_, &vars_);
	nvars_ = 0;
	naux_ = 0;
	probability_ = NULL;

	FREE_ARRAY_PTR(recvcounts_);
	FREE_ARRAY_PTR(displs_);
	FREE_ARRAY_PTR(cut_indices_);
	FREE_ARRAY_PTR(cut_status_);

	return SCIP_OKAY;
}

/** transforms constraint data into data belonging to the transformed problem */
SCIP_DECL_CONSTRANS(SCIPconshdlrBendersWorker::scip_trans)
{
	DSPdebugMessage("scip_trans\n");
	SCIP_CONSDATA * targetdata = NULL;

	//SCIP_CALL(SCIPallocMemory(scip, &targetdata));
	DSPdebugMessage("cons %p name: %s\n", sourcecons, SCIPconsGetName(sourcecons));

	/* create target constraint */
	SCIP_CALL(
			SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons),
					conshdlr, targetdata, SCIPconsIsInitial(sourcecons),
					SCIPconsIsSeparated(sourcecons),
					SCIPconsIsEnforced(sourcecons),
					SCIPconsIsChecked(sourcecons),
					SCIPconsIsPropagated(sourcecons),
					SCIPconsIsLocal(sourcecons),
					SCIPconsIsModifiable(sourcecons),
					SCIPconsIsDynamic(sourcecons),
					SCIPconsIsRemovable(sourcecons),
					SCIPconsIsStickingAtNode(sourcecons)));

	return SCIP_OKAY;
}

/** separation method of constraint handler for LP solution */
SCIP_DECL_CONSSEPALP(SCIPconshdlrBendersWorker::scip_sepalp)
{
	*result = SCIP_DIDNOTFIND;
	SCIP_CALL(sepaBenders(scip, conshdlr, NULL, from_scip_sepalp, result));
	DSPdebugMessage("scip_sepalp results in %d stage(%d)\n", *result, SCIPgetStage(scip));
	return SCIP_OKAY;
}

/** separation method of constraint handler for LP solution */
SCIP_DECL_CONSSEPASOL(SCIPconshdlrBendersWorker::scip_sepasol)
{
	*result = SCIP_DIDNOTFIND;
	SCIP_CALL(sepaBenders(scip, conshdlr, NULL, from_scip_sepasol, result));
	DSPdebugMessage("scip_sepasol results in %d stage(%d)\n", *result, SCIPgetStage(scip));
	return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for LP solutions */
SCIP_DECL_CONSENFOLP(SCIPconshdlrBendersWorker::scip_enfolp)
{
	*result = SCIP_FEASIBLE;
	SCIP_CALL(sepaBenders(scip, conshdlr, NULL, from_scip_enfolp, result));
	DSPdebugMessage("scip_enfolp results in %d stage(%d)\n", *result, SCIPgetStage(scip));
	return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for pseudo solutions */
SCIP_DECL_CONSENFOPS(SCIPconshdlrBendersWorker::scip_enfops)
{
	*result = SCIP_FEASIBLE;
	SCIP_CALL(sepaBenders(scip, conshdlr, NULL, from_scip_enfops, result));
	if (*result == SCIP_SEPARATED) *result = SCIP_INFEASIBLE;
	DSPdebugMessage("scip_enfops results in %d stage(%d)\n", *result, SCIPgetStage(scip));
	return SCIP_OKAY;
}

/** feasibility check method of constraint handler for primal solutions */
SCIP_DECL_CONSCHECK(SCIPconshdlrBendersWorker::scip_check)
{
	*result = SCIP_FEASIBLE;
	SCIP_CALL(sepaBenders(scip, conshdlr, sol, from_scip_check, result));
	DSPdebugMessage("scip_check results in %d stage(%d)\n", *result, SCIPgetStage(scip));
	return SCIP_OKAY;
}

/** variable rounding lock method of constraint handler */
SCIP_DECL_CONSLOCK(SCIPconshdlrBendersWorker::scip_lock)
{
	DSPdebugMessage("scip_lock\n");
	for (int j = 0; j < nvars_; ++j)
		SCIP_CALL(SCIPaddVarLocks(scip, vars_[j], nlockspos + nlocksneg, nlockspos + nlocksneg));

	return SCIP_OKAY;
}

/** clone method which will be used to copy constraint handler and variable pricer objects */
SCIP_DECL_CONSHDLRCLONE(scip::ObjProbCloneable* SCIPconshdlrBendersWorker::clone)
{
	*valid = true;
	SCIPconshdlrBendersWorker * conshdlrclone = new SCIPconshdlrBendersWorker(scip, scip_sepapriority_, comm_);
	conshdlrclone->setOriginalVariables(nvars_, vars_);
	return conshdlrclone;
}

SCIP_RETCODE SCIPconshdlrBendersWorker::sepaBenders(
		SCIP * scip,
		SCIP_CONSHDLR * conshdlr,
		SCIP_SOL * sol,
		whereFrom where,
		SCIP_RESULT * result)
{
	OsiCuts cuts, cuts_collected; /**< Benders cut placeholder */
	CoinPackedVector vec;
	SCIP_Real * vals = NULL; /**< current solution */

#if 0
	if (scip_checkpriority_ < 0)
	{
		/** consider incumbent solutions only */
		double primObj = SCIPgetPrimalbound(scip);
		double currObj = SCIPgetSolOrigObj(scip, sol);
		if (SCIPisLT(scip, primObj, currObj))
		{
			DSPdebugMessage(" -> primObj %e currObj %e\n", primObj, currObj);
			if (where != from_scip_check)
				*result = SCIP_DIDNOTRUN;
			return SCIP_OKAY;
		}
	}
#endif

	/** allocate memory */
	SCIP_CALL(SCIPallocMemoryArray(scip, &vals, nvars_));

	/** get current solution */
	SCIP_CALL(SCIPgetSolVals(scip, sol, nvars_, vars_, vals));

	/** TODO Does this work? */
#if 0
	if (where != from_scip_check &&
		where != from_scip_enfolp &&
		where != from_scip_enfops)
	{
		double maxviol = 1.e-10;
		for (int j = 0; j < nvars_ - naux_; ++j)
		{
			SCIP_VARTYPE vartype = SCIPvarGetType(vars_[j]);
			if (vartype == SCIP_VARTYPE_CONTINUOUS) continue;

			double viol = 0.5 - fabs(vals[j] - floor(vals[j]) - 0.5);
			//DSPdebugMessage("var[%d] %e\n", j, vals[j]);
			if (viol > maxviol)
				maxviol = viol;
		}
		DSPdebugMessage("where %d maxviol %e\n", where, maxviol);

		if (maxviol > 1.e-7)
		{
			*result = SCIP_DIDNOTRUN;
			/** free memory */
			SCIPfreeMemoryArray(scip, &vals);
			return SCIP_OKAY;
		}
	}
#endif

	if (SCIPgetStage(scip) == SCIP_STAGE_SOLVING ||
		SCIPgetStage(scip) == SCIP_STAGE_SOLVED ||
		SCIPgetStage(scip) == SCIP_STAGE_EXITSOLVE)
	{
		bool addedPoolCut = false;
		int numPoolCuts = SCIPgetNPoolCuts(scip);
		int numCutsToScan = 100;
		SCIP_CUT ** poolcuts = SCIPgetPoolCuts(scip);
		for (int i = numPoolCuts - 1; i >= 0; --i)
		{
			if (i < 0) break;
			if (numCutsToScan == 0) break;

			/** retrieve row */
			SCIP_ROW * poolcutrow = SCIPcutGetRow(poolcuts[i]);

			/** benders? */
			if (strcmp(SCIProwGetName(poolcutrow), "benders") != 0)
				continue;

			/** counter */
			numCutsToScan--;

			if (SCIPgetCutEfficacy(scip, sol, poolcutrow) > 1.e-6)
			{
				if (where == from_scip_sepalp ||
					where == from_scip_sepasol ||
					where == from_scip_enfolp)
				{
					/** add cut */
					SCIP_Bool infeasible;
					SCIP_CALL(SCIPaddCut(scip, sol, poolcutrow,
							FALSE, /**< force cut */
							&infeasible));

					if (infeasible)
						*result = SCIP_CUTOFF;
					else
						*result = SCIP_SEPARATED;
				}
				else
					*result = SCIP_INFEASIBLE;
				addedPoolCut = true;
				break;
			}
		}
		if (addedPoolCut)
		{
			DSPdebugMessage("Added pool cut\n");
			/** free memory */
			SCIPfreeMemoryArray(scip, &vals);
			return SCIP_OKAY;
		}
	}

	/** Tell workers to generate cuts */
	int message = TssBdMpi::MASTER_NEEDS_CUTS;
	MPI_Bcast(&message, 1, MPI_INT, 0, comm_);

	/** Send solutions to the workers */
	MPI_Bcast(vals, nvars_, MPI_DOUBLE, 0, comm_);

	/** Collect cut generation stutus */
	MPI_Gatherv(NULL, 0, MPI_INT, cut_status_, recvcounts_, displs_, MPI_INT, 0, comm_);

	/** Collect cuts */
	MPIgatherOsiCuts(comm_, cuts, cuts_collected);
	DSPdebugMessage("[%d]: Collected %d cuts\n", comm_rank_, cuts_collected.sizeCuts());

	/** aggregate cuts */
	aggregateCuts(cuts_collected, cuts);
	DSPdebug(cuts.printCuts());

	/** cleanup cuts_out */
	for (int i = 0; i < cuts_collected.sizeCuts(); ++i)
	{
		OsiRowCut * rc = cuts_collected.rowCutPtr(i);
		FREE_PTR(rc);
	}
	cuts_collected.dumpCuts();

	/** If found Benders cuts */
	for (int i = 0; i < cuts.sizeCuts(); ++i)
	{
		/** get cut pointer */
		OsiRowCut * rc = cuts.rowCutPtr(i);
		if (!rc) continue;

		const CoinPackedVector cutrow = rc->row();
		if (cutrow.getNumElements() == 0) continue;

		/** is optimality cut? */
		bool isOptimalityCut = false;
		for (int j = nvars_ - naux_; j < nvars_; ++j)
		{
			if (cutrow.getMaxIndex() == j)
			{
				isOptimalityCut = true;
				break;
			}
		}

		double efficacy = rc->violated(vals) / cutrow.twoNorm();
		SCIP_Bool isEfficacious = efficacy > 1.e-6;

		if (SCIPgetStage(scip) == SCIP_STAGE_INITSOLVE ||
			SCIPgetStage(scip) == SCIP_STAGE_SOLVING)
		{
			/** create empty row */
			SCIP_ROW * row = NULL;
			SCIP_CALL(SCIPcreateEmptyRowCons(scip, &row, conshdlr, "benders", rc->lb(), SCIPinfinity(scip),
					FALSE, /**< is row local? */
					FALSE, /**< is row modifiable? */
					FALSE  /**< is row removable? can this be TRUE? */));

			/** cache the row extension and only flush them if the cut gets added */
			SCIP_CALL(SCIPcacheRowExtensions(scip, row));

			/** collect all non-zero coefficients */
			for (int j = 0; j < cutrow.getNumElements(); ++j)
				SCIP_CALL(SCIPaddVarToRow(scip, row, vars_[cutrow.getIndices()[j]], cutrow.getElements()[j]));

			DSPdebugMessage("found Benders (%s) cut: act=%f, lhs=%f, norm=%f, eff=%f, min=%f, max=%f (range=%f)\n",
				isOptimalityCut ? "opti" : "feas",
				SCIPgetRowLPActivity(scip, row), SCIProwGetLhs(row), SCIProwGetNorm(row),
				SCIPgetCutEfficacy(scip, sol, row),
				SCIPgetRowMinCoef(scip, row), SCIPgetRowMaxCoef(scip, row),
				SCIPgetRowMaxCoef(scip, row)/SCIPgetRowMinCoef(scip, row));

			/** flush all changes before adding cut */
			SCIP_CALL(SCIPflushRowExtensions(scip, row));

			DSPdebugMessage("efficacy %e isEfficatious %d\n", efficacy, isEfficacious);

			if (isEfficacious)
			{
				if (where == from_scip_sepalp ||
					where == from_scip_sepasol ||
					where == from_scip_enfolp)
				{
					/** add cut */
					SCIP_Bool infeasible;
					SCIP_CALL(SCIPaddCut(scip, sol, row,
							FALSE, /**< force cut */
							&infeasible));

					if (infeasible)
						*result = SCIP_CUTOFF;
					else
						*result = SCIP_SEPARATED;
				}
				else
					*result = SCIP_INFEASIBLE;
			}

			/** add cut to global pool */
			SCIP_CALL(SCIPaddPoolCut(scip, row));
			DSPdebugMessage("number of cuts in global cut pool: %d\n", SCIPgetNPoolCuts(scip));

			/** release the row */
			SCIP_CALL(SCIPreleaseRow(scip, &row));
		}
		else if (isEfficacious &&
					where != from_scip_sepalp &&
					where != from_scip_sepasol &&
					where != from_scip_enfolp)
			*result = SCIP_INFEASIBLE;
	}

	/** free memory */
	SCIPfreeMemoryArray(scip, &vals);

	return SCIP_OKAY;
}

/** aggregate cuts in cuts_in to cuts_out */
void SCIPconshdlrBendersWorker::aggregateCuts(OsiCuts cuts_in, OsiCuts &cuts_out)
{
#define FREE_MEMORY                  \
	FREE_2D_ARRAY_PTR(naux_,aggval); \
	FREE_ARRAY_PTR(aggrhs);

	double ** aggval  = NULL; /** aggregated dense cut coefficients */
	double *  aggrhs  = NULL; /** aggregated cut rhs */
	CoinPackedVector vec;
	bool isInfeasible = false;

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

	/** cleanup cuts_out */
	for (int i = 0; i < cuts_out.sizeCuts(); ++i)
	{
		OsiRowCut * rc = cuts_out.rowCutPtr(i);
		FREE_PTR(rc);
	}
	cuts_out.dumpCuts();

	/** aggregate cuts */
	for (int i = 0; i < cuts_in.sizeCuts(); ++i)
	{
		int s = cut_indices_[i];
		double probability = probability_ ? probability_[s] : 1.0;
		OsiRowCut * rc = cuts_in.rowCutPtr(i);

		/** feasibility cut? */
		DSPdebugMessage("cut_indices %d cut_status %d weight %e\n", s, cut_status_[s], probability);
		if (cut_status_[s] == STO_STAT_PRIM_INFEASIBLE)
		{
			/** initialize vector */
			vec.clear();

			/** set row vector */
			vec = probability * rc->row();

			OsiRowCut fcut;
			fcut.setRow(vec);
			fcut.setUb(COIN_DBL_MAX);
			fcut.setLb(probability * rc->lb());

			cuts_out.insert(fcut);
			isInfeasible = true;
			break;
		}

		/** consider only optimality cuts */
		if (cut_status_[s] != STO_STAT_OPTIMAL)
			continue;

		/** calculate weighted aggregation of cuts */
		const CoinPackedVector row = rc->row();
		int ind_aux = s % naux_;
		for (int j = 0; j < row.getNumElements(); ++j)
			aggval[ind_aux][row.getIndices()[j]] += probability * row.getElements()[j];
		aggrhs[ind_aux] += probability * rc->lb();
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
			{
				if (fabs(aggval[s][j]) > 1E-10)
					vec.insert(j, aggval[s][j]);
			}

			if (fabs(aggrhs[s]) < 1E-10)
				aggrhs[s] = 0.0;

			OsiRowCut rc;
			rc.setRow(vec);
			rc.setUb(COIN_DBL_MAX); /** TODO: for minimization */
			rc.setLb(aggrhs[s]);

			cuts_out.insert(rc);
		}
	}

	END_TRY_CATCH(FREE_MEMORY)

	FREE_MEMORY

#undef FREE_MEMORY
}

void SCIPconshdlrBendersWorker::init(
		int nsubprobs,             /**< number of subproblems */
		const double * probability /**< probability of each subproblem */)
{
	recvcounts_  = new int [comm_size_];
	displs_      = new int [comm_size_];
	probability_ = probability;
	cut_indices_ = new int [nsubprobs];
	cut_status_  = new int [nsubprobs];

	/** The master does not solve subproblem. */
	recvcounts_[0] = 0;
	displs_[0] = 0;

	/** Subproblems are assigned to each process in round-and-robin fashion. */
	for (int i = 1, j = 0; i < comm_size_; ++i)
	{
		recvcounts_[i] = 0;
		for (int s = i; s < nsubprobs; s += comm_size_-1)
		{
			recvcounts_[i]++;
			cut_indices_[j++] = s;
		}
		displs_[i] = i == 0 ? 0 : displs_[i-1] + recvcounts_[i-1];
		DSPdebugMessage("recvcounts_[%d] %d displs_[%d] %d\n", i, recvcounts_[i], i, displs_[i]);
	}
}
