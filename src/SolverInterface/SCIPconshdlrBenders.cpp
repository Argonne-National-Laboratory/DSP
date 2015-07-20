/*
 * SCIPconshdlrBenders.cpp
 *
 *  Created on: Dec 8, 2014
 *      Author: kibaekkim
 */

//#define DSP_DEBUG

/** DSP */
#include "SolverInterface/SCIPconshdlrBenders.h"
#include "Utility/StoMessage.h"

/** constraint data for Benders cuts */
struct SCIP_ConsData
{
	/** nothing */
};

/** destructor of constraint handler to free user data (called when SCIP is exiting) */
SCIP_DECL_CONSFREE(SCIPconshdlrBenders::scip_free)
{
	for (int j = 0; j < nvars_; ++j)
		vars_[j] = NULL;
	SCIPfreeMemoryArray(scip_, &vars_);
	nvars_ = 0;
	naux_ = 0;

	return SCIP_OKAY;
}

/** transforms constraint data into data belonging to the transformed problem */
SCIP_DECL_CONSTRANS(SCIPconshdlrBenders::scip_trans)
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
SCIP_DECL_CONSSEPALP(SCIPconshdlrBenders::scip_sepalp)
{
	*result = SCIP_DIDNOTFIND;
	SCIP_CALL(sepaBenders(scip, conshdlr, NULL, from_scip_sepalp, result));
	DSPdebugMessage("scip_sepalp results in %d\n", *result);
	return SCIP_OKAY;
}

/** separation method of constraint handler for LP solution */
SCIP_DECL_CONSSEPASOL(SCIPconshdlrBenders::scip_sepasol)
{
	*result = SCIP_DIDNOTFIND;
	SCIP_CALL(sepaBenders(scip, conshdlr, NULL, from_scip_sepasol, result));
	DSPdebugMessage("scip_sepasol results in %d\n", *result);
	return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for LP solutions */
SCIP_DECL_CONSENFOLP(SCIPconshdlrBenders::scip_enfolp)
{
	*result = SCIP_FEASIBLE;
	SCIP_CALL(sepaBenders(scip, conshdlr, NULL, from_scip_enfolp, result));
	DSPdebugMessage("scip_enfolp results in %d\n", *result);
	return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for pseudo solutions */
SCIP_DECL_CONSENFOPS(SCIPconshdlrBenders::scip_enfops)
{
	*result = SCIP_FEASIBLE;
	SCIP_CALL(sepaBenders(scip, conshdlr, NULL, from_scip_enfops, result));
	if (*result == SCIP_SEPARATED) *result = SCIP_INFEASIBLE;
	DSPdebugMessage("scip_enfops results in %d\n", *result);
	return SCIP_OKAY;
}

/** feasibility check method of constraint handler for primal solutions */
SCIP_DECL_CONSCHECK(SCIPconshdlrBenders::scip_check)
{
	*result = SCIP_FEASIBLE;
	SCIP_CALL(sepaBenders(scip, conshdlr, sol, from_scip_check, result));
	DSPdebugMessage("scip_check results in %d\n", *result);
	return SCIP_OKAY;
}

/** variable rounding lock method of constraint handler */
SCIP_DECL_CONSLOCK(SCIPconshdlrBenders::scip_lock)
{
	DSPdebugMessage("scip_lock\n");
	for (int j = 0; j < nvars_; ++j)
		SCIP_CALL(SCIPaddVarLocks(scip, vars_[j], nlockspos + nlocksneg, nlockspos + nlocksneg));

	return SCIP_OKAY;
}

SCIP_RETCODE SCIPconshdlrBenders::sepaBenders(
		SCIP * scip,
		SCIP_CONSHDLR * conshdlr,
		SCIP_SOL * sol,
		whereFrom where,
		SCIP_RESULT * result)
{
	OsiCuts cs; /**< Benders cut placeholder */
	SCIP_Real * vals = NULL; /**< current solution */

#if 1
	if (scip_checkpriority_ < 0)
	{
		/** consider incumbent solutions only */
		double primObj = SCIPgetPrimalbound(scip);
		double currObj = SCIPgetSolOrigObj(scip, sol);
		if (SCIPisLT(scip, primObj, currObj))
		{
			DSPdebugMessage(" -> primObj %e currObj %e\n", primObj, currObj);
			return SCIP_OKAY;
		}
	}
#endif

	/** allocate memory */
	SCIP_CALL(SCIPallocMemoryArray(scip, &vals, nvars_));

	/** get current solution */
	SCIP_CALL(SCIPgetSolVals(scip, sol, nvars_, vars_, vals));

	/** TODO The following filter does not work, meaning that it provides suboptimal solution.
	 * I do not know the reason. */
#if 0
	double maxviol = 1.e-10;
	for (int j = 0; j < nvars_ - naux_; ++j)
	{
		SCIP_VARTYPE vartype = SCIPvarGetType(vars_[j]);
		if (vartype == SCIP_VARTYPE_CONTINUOUS) continue;

		double viol = 0.5 - fabs(vals[j] - floor(vals[j]) - 0.5);
		if (viol > maxviol)
			maxviol = viol;
	}
	DSPdebugMessage("maximum violation %e\n", maxviol);

	if (where != from_scip_check &&
		where != from_scip_enfolp &&
		where != from_scip_enfops &&
		maxviol > 1.e-7)
	{
		printf("where %d maxviol %e\n", where, maxviol);
		/** free memory */
		SCIPfreeMemoryArray(scip, &vals);
		return SCIP_OKAY;
	}
#endif

#ifdef DSP_DEBUG2
	double minvals = COIN_DBL_MAX;
	double maxvals = -COIN_DBL_MAX;
	double sumvals = 0.;
	double ssvals  = 0.;
	//printf("nvars_ %d naux_ %d nAuxvars_ %d\n", nvars_, naux_, tss_->nAuxvars_);
	for (int j = 0; j < nvars_ - naux_; ++j)
	{
//		if (vals[j] < 0 || vals[j] > 1)
//			printf("solution %d has value %e.\n", j, vals[j]);
		sumvals += vals[j];
		ssvals  += vals[j] * vals[j];
		minvals = minvals > vals[j] ? vals[j] : minvals;
		maxvals = maxvals < vals[j] ? vals[j] : maxvals;
	}
	DSPdebugMessage("solution: min %e max %e avg %e sum %e two-norm %e\n",
			minvals, maxvals, sumvals / nvars_, sumvals, sqrt(ssvals));
#endif

#define SCAN_GLOBAL_CUT_POOL
#ifdef SCAN_GLOBAL_CUT_POOL
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
					else //if (*result != SCIP_CUTOFF)
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
#endif

	/** generate Benders cuts */
	assert(tss_);
	tss_->generateCuts(nvars_, vals, &cs);

	/** If found Benders cuts */
	for (int i = 0; i < cs.sizeCuts(); ++i)
	{
		/** get cut pointer */
		OsiRowCut * rc = cs.rowCutPtr(i);
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

#define KK_TEST
#ifdef KK_TEST
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
					else //if (*result != SCIP_CUTOFF)
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
#else
		if (where == from_scip_sepalp ||
			where == from_scip_sepasol ||
			where == from_scip_enfolp)
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
				SCIPgetCutEfficacy(scip, NULL, row),
				SCIPgetRowMinCoef(scip, row), SCIPgetRowMaxCoef(scip, row),
				SCIPgetRowMaxCoef(scip, row)/SCIPgetRowMinCoef(scip, row));

			/** flush all changes before adding cut */
			SCIP_CALL(SCIPflushRowExtensions(scip, row));

			/** is cut efficacious? */
			if (isOptimalityCut)
			{
				efficacy = SCIPgetCutEfficacy(scip, sol, row);
				isEfficacious = SCIPisCutEfficacious(scip, sol, row);
			}
			else
			{
				efficacy = rc->violated(vals);
				isEfficacious = efficacy > 1.e-6;
			}

			if (isEfficacious)
			{
				/** add cut */
				SCIP_Bool infeasible;
				SCIP_CALL(SCIPaddCut(scip, sol, row,
						FALSE, /**< force cut */
						&infeasible));

				if (infeasible)
					*result = SCIP_CUTOFF;
				else if (*result != SCIP_CUTOFF)
					*result = SCIP_SEPARATED;
			}

			/** add cut to global pool */
			SCIP_CALL(SCIPaddPoolCut(scip, row));

			/** release the row */
			SCIP_CALL(SCIPreleaseRow(scip, &row));
		}
		else
		{
			if (isOptimalityCut)
			{
				efficacy = rc->violated(vals) / cutrow.twoNorm();
				isEfficacious = efficacy > 0.05;
			}
			else
			{
				efficacy = rc->violated(vals);
				isEfficacious = efficacy > 1.e-6;
			}
			DSPdebugMessage("%s efficacy %e\n", isOptimalityCut ? "Opti" : "Feas", efficacy);

			if (isEfficacious == TRUE)
				*result = SCIP_INFEASIBLE;
		}
#endif
	}

	/** free memory */
	SCIPfreeMemoryArray(scip, &vals);

	return SCIP_OKAY;
}

/** set original variable pointers */
SCIP_RETCODE SCIPconshdlrBenders::setOriginalVariables(int nvars, SCIP_Var ** vars)
{
	nvars_ = nvars;
	SCIP_CALL(SCIPallocMemoryArray(scip_, &vars_, nvars_));
	for (int j = 0; j < nvars_; ++j)
		vars_[j] = vars[j];
	return SCIP_OKAY;
}

/** clone method which will be used to copy constraint handler and variable pricer objects */
SCIP_DECL_CONSHDLRCLONE(scip::ObjProbCloneable* SCIPconshdlrBenders::clone)
{
	*valid = true;
#if 1
	SCIPconshdlrBenders * conshdlrclone = new SCIPconshdlrBenders(scip, scip_sepapriority_);
	conshdlrclone->assignTssBdSub(tss_);
	conshdlrclone->setOriginalVariables(nvars_, vars_);
	return conshdlrclone;
#else
	return NULL;
#endif
}

/** creates and captures a Benders constraint */
SCIP_RETCODE SCIPcreateConsBenders(
		SCIP * scip,
		SCIP_CONS ** cons,
		const char * name)
{

	SCIP_CONSHDLR * conshdlr = NULL;
	SCIP_CONSDATA * consdata = NULL;

	/* find the Benders constraint handler */
	conshdlr = SCIPfindConshdlr(scip, name);
	assert(conshdlr);

	/* create constraint data */
	//SCIP_CALL(SCIPallocMemory(scip_, &consdata));

	/* create constraint */
	SCIP_CALL(SCIPcreateCons(scip, cons, name, conshdlr, consdata,
			false, /**< being in the initial LP? */
			true,  /**< separated during LP process? */
			true,  /**< enforced? */
			true,  /**< checked for feasibility? */
			true,  /**< propagate? */
			false, /**< only locally valid? */
			false, /**< modifiable? */
			false, /**< is constraint subject to aging? */
			true,  /**< removable? */
			false));

	return SCIP_OKAY;
}

