/*
 * SCIPconshdlrBendersDd.cpp
 *
 *  Created on: Dec 14, 2014
 *      Author: kibaekkim
 */

//#define DSP_DEBUG

#include "SolverInterface/SCIPconshdlrBendersDd.h"

#include <Utility/DspMessage.h>

#define DSP_COLLECT_CUTS 1
#define KK_DEBUG_PRINT 0

/** constraint data for Benders cuts */
struct SCIP_ConsData
{
	/** nothing to do */
};

/** destructor of constraint handler to free user data (called when SCIP is exiting) */
SCIP_DECL_CONSFREE(SCIPconshdlrBendersDd::scip_free)
{
	DSPdebugMessage("scip_free\n");

	/** release variables */
	for (int j = 0; j < nvars_; ++j)
		vars_[j] = NULL;
	SCIPfreeMemoryArray(scip_, &vars_);

	/** release cut pool */
#if 0
	for (int i = 0; i < cutsToAdd_->sizeCuts(); ++i)
	{
		OsiRowCut * rc = cutsToAdd_->rowCutPtr(i);
		FREE_PTR(rc);
	}
#endif
	cutsToAdd_->dumpCuts();
	for (int i = 0; i < cutsAdded_->sizeCuts(); ++i)
	{
		OsiRowCut * rc = cutsAdded_->rowCutPtr(i);
		FREE_PTR(rc);
	}
	cutsAdded_->dumpCuts();

	nvars_ = 0;
	ncols_ = 0;
	ncols_first_ = 0;
	obj_first_ = NULL;

	return SCIP_OKAY;
}

/** transforms constraint data into data belonging to the transformed problem */
SCIP_DECL_CONSTRANS(SCIPconshdlrBendersDd::scip_trans)
{
	DSPdebugMessage("scip_trans\n");
	SCIP_CONSDATA * targetdata = NULL;

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
SCIP_DECL_CONSSEPALP(SCIPconshdlrBendersDd::scip_sepalp)
{
	*result = SCIP_DIDNOTFIND;
	SCIP_CALL(sepaBenders(scip, conshdlr, NULL, from_scip_sepalp, result));
	DSPdebugMessage("scip_sepalp results in %d\n", *result);
	return SCIP_OKAY;
}

/** separation method of constraint handler for LP solution */
SCIP_DECL_CONSSEPASOL(SCIPconshdlrBendersDd::scip_sepasol)
{
	*result = SCIP_DIDNOTFIND;
	SCIP_CALL(sepaBenders(scip, conshdlr, NULL, from_scip_sepasol, result));
	DSPdebugMessage("scip_sepasol results in %d\n", *result);
	return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for LP solutions */
SCIP_DECL_CONSENFOLP(SCIPconshdlrBendersDd::scip_enfolp)
{
	*result = SCIP_FEASIBLE;
	SCIP_CALL(sepaBenders(scip, conshdlr, NULL, from_scip_enfolp, result));
	DSPdebugMessage("%p scip_enfolp results in %d\n", scip, *result);
	return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for pseudo solutions */
SCIP_DECL_CONSENFOPS(SCIPconshdlrBendersDd::scip_enfops)
{
	*result = SCIP_FEASIBLE;
	SCIP_CALL(sepaBenders(scip, conshdlr, NULL, from_scip_enfops, result));
	if (*result == SCIP_SEPARATED) *result = SCIP_SOLVELP;
	DSPdebugMessage("scip_enfops results in %d\n", *result);
	return SCIP_OKAY;
}

/** feasibility check method of constraint handler for primal solutions */
SCIP_DECL_CONSCHECK(SCIPconshdlrBendersDd::scip_check)
{
	*result = SCIP_FEASIBLE;
	SCIP_CALL(sepaBenders(scip, conshdlr, sol, from_scip_check, result));
	DSPdebugMessage("scip_check results in %d\n", *result);
	return SCIP_OKAY;
}

/** variable rounding lock method of constraint handler */
SCIP_DECL_CONSLOCK(SCIPconshdlrBendersDd::scip_lock)
{
	DSPdebugMessage("scip_lock\n");
	for (int j = 0; j < nvars_; ++j)
		SCIP_CALL(SCIPaddVarLocks(scip, vars_[j], nlockspos + nlocksneg, nlockspos + nlocksneg));

	return SCIP_OKAY;
}

SCIP_RETCODE SCIPconshdlrBendersDd::sepaBenders(
		SCIP * scip,
		SCIP_CONSHDLR * conshdlr,
		SCIP_SOL * sol,
		whereFrom where,
		SCIP_RESULT * result)
{
	OsiCuts cs; /**< Benders cut placeholder */
	SCIP_Real * vals = NULL; /**< current solution */

#if 0
	/** TODO Consider the root node only */
	DSPdebugMessage("Node depth %d\n", SCIPgetDepth(scip));
	if (SCIPgetDepth(scip) > 0)
		return SCIP_OKAY;
#endif

	/** consider incumbent solutions only */
	double primObj = SCIPgetPrimalbound(scip);
	double currObj = SCIPgetSolOrigObj(scip, sol);
	if (SCIPisLT(scip, primObj, currObj))
	{
		DSPdebugMessage("primObj %e currObj %e\n", primObj, currObj);
		return SCIP_OKAY;
	}

	/** allocate memory */
	SCIP_CALL(SCIPallocMemoryArray(scip, &vals, nvars_));

	/** get current solution */
	SCIP_CALL(SCIPgetSolVals(scip, sol, nvars_, vars_, vals));

#define SCAN_GLOBAL_CUT_POOL
#ifdef SCAN_GLOBAL_CUT_POOL
	if (SCIPgetStage(scip) == SCIP_STAGE_SOLVING ||
		SCIPgetStage(scip) == SCIP_STAGE_SOLVED ||
		SCIPgetStage(scip) == SCIP_STAGE_EXITSOLVE)
	{
		bool addedPoolCut = false;
		int numPoolCuts = SCIPgetNPoolCuts(scip);
		int numCutsToScan = tss_ == NULL ? numPoolCuts : 100;
		SCIP_CUT ** poolcuts = SCIPgetPoolCuts(scip);
		for (int i = numPoolCuts - 1; i >= 0; --i)
		{
			if (i < 0) break;
			if (numCutsToScan == 0) break;

			/** retrieve row */
			SCIP_ROW * poolcutrow = SCIPcutGetRow(poolcuts[i]);

			/** benders? */
			if (strcmp(SCIProwGetName(poolcutrow), "bendersDd") != 0)
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

	/** update cut efficacy in cutsToAdd_ */
	bool isCutFromPool = false;
	int maxEfficaciousCut = -1;
	double maxEfficacy = 1.e-6;
	DSPdebugMessage("Number of cutsToAdd_ %d\n", cutsToAdd_->sizeCuts());
	for (int i = cutsToAdd_->sizeCuts() - 1; i >= 0; --i)
	{
		OsiRowCut * rc = cutsToAdd_->rowCutPtr(i);
		if (!rc) continue;
		const CoinPackedVector row = rc->row();

		/** is optimality cut? */
		DSPdebugMessage("row.getNumElements() %d\n", row.getNumElements());
		bool isOptimalityCut = row.getIndices()[row.getNumElements() - 1] == nvars_ - 1;

		/** calculate efficacy */
		double efficacy = rc->violated(vals);
		if (isOptimalityCut) efficacy /= row.twoNorm();

		/** determine if efficacious */
		if (efficacy <= 1.e-6)
			continue;
		else if (efficacy > maxEfficacy)
		{
			maxEfficaciousCut = i;
			maxEfficacy = efficacy;
		}
	}

	DSPdebugMessage("maxEfficacy %e\n", maxEfficacy);
	if (maxEfficaciousCut >= 0)
	{
		/** move one from cut pool */
		OsiRowCut * rc = NULL;
		if (where == from_scip_check)
			rc = cutsToAdd_->rowCutPtr(maxEfficaciousCut);
		else
			rc = cutsToAdd_->rowCutPtrAndZap(maxEfficaciousCut);
		cs.insert(rc);

		/** mark as cut from pool */
		isCutFromPool = true;

		/** increment counter */
		nCutsReused_++;
	}
	else if (tss_ != NULL && SCIPgetDepth(scip) == 0)
	{
		/** generate Benders cuts */
		tss_->generateCuts(nvars_, vals, &cs, TssBdSub::TssDd);
		DSPdebug(cs.printCuts());

		/** construct upper bounding cuts */
		constructCuts(cs);
	}

	/** if any cut to add? */
	for (int i = 0; i < cs.sizeCuts(); ++i)
	{
		/** get cut pointer */
		OsiRowCut * rc = cs.rowCutPtr(i);
		if (!rc) continue;

		const CoinPackedVector cutrow = rc->row();
		if (cutrow.getNumElements() == 0) continue;

		/** is optimality cut? */
		bool isOptimalityCut = cutrow.getIndices()[cutrow.getNumElements() - 1] == nvars_ - 1;

		/** is cut efficacious? */
		double efficacy;
		SCIP_Bool isEfficacious;
		if (isCutFromPool)
		{
			efficacy = maxEfficacy;
			isEfficacious = true;
		}
		else
		{
			efficacy = rc->violated(vals) / cutrow.twoNorm();
			isEfficacious = efficacy > 1.e-6;
		}
		DSPdebugMessage("efficacy %e\n", efficacy);

		if (SCIPgetStage(scip) == SCIP_STAGE_INITSOLVE ||
			SCIPgetStage(scip) == SCIP_STAGE_SOLVING)
		{
			/** create empty row */
			SCIP_ROW * row = NULL;
			SCIP_CALL(SCIPcreateEmptyRowCons(scip, &row, conshdlr, "bendersDd", rc->lb(), SCIPinfinity(scip),
					FALSE, /**< is row local? */
					FALSE, /**< is row modifiable? */
					FALSE  /**< is row removable? can this be TRUE? */));

			/** cache the row extension and only flush them if the cut gets added */
			SCIP_CALL(SCIPcacheRowExtensions(scip, row));

			/** collect all non-zero coefficients */
			for (int j = 0; j < cutrow.getNumElements(); ++j)
				SCIP_CALL(SCIPaddVarToRow(scip, row, vars_[cutrow.getIndices()[j]], cutrow.getElements()[j]));

			/** flush all changes before adding cut */
			SCIP_CALL(SCIPflushRowExtensions(scip, row));

			DSPdebugMessage("found Benders cut (%s): act=%f, lhs=%f, norm=%f, eff=%f (%s), min=%f, max=%f (range=%f)\n",
					isOptimalityCut ? "opti" : "feas",
					sol ? SCIPgetRowSolActivity(scip, row, sol) : SCIPgetRowLPActivity(scip, row),
					SCIProwGetLhs(row), SCIProwGetNorm(row),
					efficacy, isEfficacious ? "yes" : "no",
					SCIPgetRowMinCoef(scip, row), SCIPgetRowMaxCoef(scip, row),
					SCIPgetRowMaxCoef(scip, row)/SCIPgetRowMinCoef(scip, row));

			if (where == from_scip_sepalp ||
				where == from_scip_sepasol ||
				where == from_scip_enfolp)
			{
				if (isEfficacious)
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

				if (!isCutFromPool)
					cutsAdded_->insert(rc);
			}
			else
			{
				if (isEfficacious)
					*result = SCIP_INFEASIBLE;
				if (!isCutFromPool)
					FREE_PTR(rc);
			}

			/** add cut to global pool */
			if (SCIPgetStage(scip) == SCIP_STAGE_SOLVING)
				SCIP_CALL(SCIPaddPoolCut(scip, row));

			/** release the row */
			SCIP_CALL(SCIPreleaseRow(scip, &row));
		}
		else
		{
			if (where == from_scip_sepalp ||
				where == from_scip_sepasol ||
				where == from_scip_enfolp)
			{
				if (!isCutFromPool)
					cutsAdded_->insert(rc);
			}
			else
			{
				if (isEfficacious)
					*result = SCIP_INFEASIBLE;
				if (!isCutFromPool)
					FREE_PTR(rc);
			}
		}
	}

	/** dump local cutpool */
	cs.dumpCuts();

	/** free memory */
	SCIPfreeMemoryArray(scip, &vals);

	return SCIP_OKAY;
}

/** construct upper bounding cuts */
void SCIPconshdlrBendersDd::constructCuts(OsiCuts & cs)
{
	assert(obj_first_);
	OsiCuts cuts;

	for (int i = 0; i < cs.sizeCuts(); ++i)
	{
		/** retrieve cut data */
		OsiRowCut *      rc = cs.rowCutPtr(i);
		if (rc == NULL) continue;
		//rc->print();

		CoinPackedVector row = rc->row();
		double *         dense  = row.denseVector(ncols_);
		double           cutrhs = rc->lb();
		bool isFcut = fabs(dense[ncols_ - 1]) > 1E-10 ? false : true;

		/** construct upper bounding cuts */
		row.clear();
		for (int j = 0; j < ncols_first_; ++j)
		{
			if (!isFcut)
				dense[j] -= obj_first_[j];
			if (fabs(dense[j]) > 1E-10)
				row.insert(j, dense[j]);
		}
		if (!isFcut)
			row.insert(ncols_ - 1, 1.0);

		/** create new bounding cut */
		OsiRowCut rowcut;
		rowcut.setRow(row);
		rowcut.setUb(COIN_DBL_MAX);
		rowcut.setLb(cutrhs);
		//rowcut.print();

		/** store cut */
		cuts.insert(rowcut);

		/** free memory */
		FREE_ARRAY_PTR(dense);
		FREE_PTR(rc);
	}

	/** clean up original cut pool */
	cs.dumpCuts();

	/** copy new cuts */
	for (int i = 0; i < cuts.sizeCuts(); ++i)
		cs.insert(cuts.rowCut(i));
}

/** set original variable pointers */
SCIP_RETCODE SCIPconshdlrBendersDd::setOriginalVariables(int nvars, SCIP_Var ** vars)
{
	nvars_ = nvars;
	SCIP_CALL(SCIPallocMemoryArray(scip_, &vars_, nvars_));
	for (int j = 0; j < nvars_; ++j)
		vars_[j] = vars[j];
	return SCIP_OKAY;
}

/** clone method which will be used to copy constraint handler and variable pricer objects */
SCIP_DECL_CONSHDLRCLONE(scip::ObjProbCloneable* SCIPconshdlrBendersDd::clone)
{
	*valid = true;
	SCIPconshdlrBendersDd * conshdlrclone = new SCIPconshdlrBendersDd(scip);
	conshdlrclone->assignData(tss_, ncols_, ncols_first_, obj_first_);
	conshdlrclone->setOriginalVariables(nvars_, vars_);
	return conshdlrclone;
}

/** set cutsToAdd_ */
void SCIPconshdlrBendersDd::setCutsToAdd(OsiCuts * cuts)
{
#if 0
	for (int i = 0; i < cutsToAdd_->sizeCuts(); ++i)
	{
		OsiRowCut * rc = cutsToAdd_->rowCutPtr(i);
		FREE_PTR(rc);
	}
	cutsToAdd_->dumpCuts();

	for (int i = 0; i < cuts->sizeCuts(); ++i)
	{
		OsiRowCut * rc = cuts->rowCutPtr(i);
		if (!rc) continue;

		/** copy */
		cutsToAdd_->insert(*rc);
	}
#else
	cutsToAdd_->dumpCuts();
	for (int i = 0; i < cuts->sizeCuts(); ++i)
	{
		OsiRowCut * rc = cuts->rowCutPtr(i);
		if (!rc) continue;

		/** shallow copy */
		cutsToAdd_->insert(rc);
	}
#endif
}

/** clean cutsAdded */
void SCIPconshdlrBendersDd::clearCutsAdded()
{
	for (int i = 0; i < cutsAdded_->sizeCuts(); ++i)
	{
		OsiRowCut * rc = cutsAdded_->rowCutPtr(i);
		FREE_PTR(rc);
	}
	cutsAdded_->dumpCuts();
}

/** creates and captures a Benders constraint */
SCIP_RETCODE SCIPcreateConsBendersDd(
		SCIP * scip,
		SCIP_CONS ** cons,
		const char * name)
{
	DSPdebugMessage("SCIPcreateConsBendersDd\n");
	SCIP_CONSHDLR * conshdlr = NULL;
	SCIP_CONSDATA * consdata = NULL;

	/* find the subtour constraint handler */
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

