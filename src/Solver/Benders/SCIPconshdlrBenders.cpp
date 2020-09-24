/*
 * SCIPconshdlrBenders.cpp
 *
 *  Created on: Dec 8, 2014
 *      Author: kibaekkim
 */

// #define DSP_DEBUG

#include "OsiCuts.hpp"

/** DSP */
#include "Utility/DspMessage.h"
#include "Model/TssModel.h"
#include "Solver/Benders/SCIPconshdlrBenders.h"

#define BENDERS_CUT_REMOVABLE TRUE
#define BENDERS_CUT_FORCE FALSE
#define INTOPT_CUT_REMOVABLE TRUE
#define INTOPT_CUT_FORCE FALSE
#define NOGOOD_CUT_REMOVABLE TRUE
#define NOGOOD_CUT_FORCE FALSE

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
	FREE_ARRAY_PTR(probability_);

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

/** constraint enforcing method of constraint handler for LP solutions */
SCIP_DECL_CONSENFOLP(SCIPconshdlrBenders::scip_enfolp)
{
#define FREE_MEMORY                 \
	FREE_ARRAY_PTR(recourse_values)

	double* recourse_values = NULL;
	double weighted_sum_of_recourse = 0.0; // weighted sum of the recourse values
	double approx_recourse = 0.0;
	double weighted_recourse;
	int nsubs = 0;
	bool is_pure_binary = (nvars_ - naux_ == SCIPgetNOrigBinVars(scip));

	/**
	 * This part computes the exact recourse function value at the current solution,
	 * 	adjusts the auxiliary variable values of the solution,
	 * 	and send it to SCIP as a primal solution.
	*/
	if (SCIPgetStage(scip) == SCIP_STAGE_SOLVING &&
			(isIntegralRecourse() || model_->isDro())) {
		SCIP_Bool stored;
		SCIP_Sol* sol;
		nsubs = bdsub_->getNumSubprobs();

		recourse_values = new double [nsubs];

		// get the current solution
		SCIP_CALL(SCIPcreateCurrentSol(scip, &sol, NULL));

		// evaluate the recourse value
		SCIP_CALL(evaluateRecourse(scip, sol, recourse_values));

		// compute weighted sum for DRO; otherwise, returns the current recourse
		computeProbability(recourse_values);
		for (int j = 0; j < nsubs; ++j) {
			DSPdebugMessage("----- scip_enfolp: exact recourse value [%d] %e with probability %e\n", j, recourse_values[j], probability_[j]);
			weighted_sum_of_recourse += recourse_values[j] * probability_[j];
		}

		// set a feasible value of the auxiliary variable
		for (int j = 0; j < naux_; ++j) {
			approx_recourse += SCIPgetSolVal(scip, sol, vars_[nvars_ - naux_ + j]);
			weighted_recourse = 0.0;
			for (int k = j; k < nsubs; k += naux_)
				weighted_recourse += recourse_values[k] * probability_[k];
			DSPdebugMessage("----- scip_enfolp: set variable[%d] from %e to %e (bounds [%e,%e])\n",
				nvars_ - naux_ + j, SCIPgetSolVal(scip, sol, vars_[nvars_ - naux_ + j]), weighted_recourse,
				SCIPvarGetLbLocal(vars_[nvars_ - naux_ + j]), 
				SCIPvarGetUbLocal(vars_[nvars_ - naux_ + j]));
			SCIP_CALL(SCIPsetSolVal(scip, sol, vars_[nvars_ - naux_ + j], weighted_recourse));
		}

		// set the solution as a primal feasible solution
		double sense = SCIPgetObjsense(scip) == SCIP_OBJSENSE_MINIMIZE ? 1.0 : -1.0;
		if (SCIPisLT(scip, (SCIPsolGetOrigObj(sol) - SCIPgetPrimalbound(scip)) * sense, 0.0)) {
			SCIP_CALL(SCIPtrySol(scip, sol, FALSE, TRUE, TRUE, TRUE, TRUE, &stored));
			DSPdebugMessage("----- scip_enfolp: is solution (obj: %e) stored? %s\n", SCIPsolGetOrigObj(sol), stored ? "yes" : "no");
		}
		SCIP_CALL(SCIPfreeSol(scip, &sol));
	}

	*result = SCIP_FEASIBLE;
	/**
	 * TODO: DRO needs to pass probability_ to sepaBenders
	 */
	SCIP_CALL(sepaBenders(scip, conshdlr, NULL, from_scip_enfolp, result));
	DSPdebugMessage("scip_enfolp results in %d stage %d\n", *result, SCIPgetStage(scip));

	if (is_pure_binary) {
		if (model_->isStochastic() && 
				isIntegralRecourse() && 
				SCIPgetStage(scip) == SCIP_STAGE_SOLVING &&
				approx_recourse < weighted_sum_of_recourse) {
			/** integer Benders cut */
			double recourse_lb = SCIPgetDualbound(scip);
			for (int j = 0; j < nvars_-naux_; ++j) {
				recourse_lb -= CoinMax(0.0, SCIPvarGetObj(vars_[j]));
			}
			addIntOptimalityCut(scip, conshdlr, weighted_sum_of_recourse, recourse_lb, result);
			DSPdebugMessage("----- scip_enfolp: integer optimality cut result %d\n", *result);
		}

		// Add no-good cut for pure binary
		addNoGoodCut(scip, conshdlr, result);
		DSPdebugMessage("----- scip_enfolp: No-good cut result %d\n", *result);
	}

	FREE_MEMORY;

	return SCIP_OKAY;
#undef FREE_MEMORY
}

/** constraint enforcing method of constraint handler for pseudo solutions */
SCIP_DECL_CONSENFOPS(SCIPconshdlrBenders::scip_enfops)
{
	*result = SCIP_FEASIBLE;
	SCIP_CALL(sepaBenders(scip, conshdlr, NULL, from_scip_enfops, result));
	if (*result == SCIP_SEPARATED) *result = SCIP_INFEASIBLE;
	return SCIP_OKAY;
}

/** feasibility check method of constraint handler for primal solutions */
SCIP_DECL_CONSCHECK(SCIPconshdlrBenders::scip_check)
{
	*result = SCIP_FEASIBLE;
	SCIP_CALL(sepaBenders(scip, conshdlr, sol, from_scip_check, result));
	DSPdebugMessage("scip_check results in %d stage(%d)\n", *result, SCIPgetStage(scip));

	if (isIntegralRecourse() && 
		SCIPgetStage(scip) == SCIP_STAGE_SOLVING &&
		*result == SCIP_FEASIBLE)
	{
		int nsubs = bdsub_->getNumSubprobs();
		double weighted_sum_of_recourse = 0.0;
		double approx_recourse = 0.0;

		double* recourse_values = new double [nsubs];

		// evaluate the recourse value
		SCIP_CALL(evaluateRecourse(scip, sol, recourse_values));

		// compute weighted sum for DRO; otherwise, returns the current recourse
		computeProbability(recourse_values);
		for (int j = 0; j < nsubs; ++j) {
			DSPdebugMessage("----- scip_check: exact recourse value [%d] %e with probability %e\n", j, recourse_values[j], probability_[j]);
			weighted_sum_of_recourse += recourse_values[j] * probability_[j];
		}

		// set a feasible value of the auxiliary variable
		for (int j = 0; j < naux_; ++j) {
			approx_recourse += SCIPgetSolVal(scip, sol, vars_[nvars_ - naux_ + j]);
		}

		if (approx_recourse < weighted_sum_of_recourse) {
			*result = SCIP_INFEASIBLE;
			DSPdebugMessage("----- scip_check: rejects solution (approx %e, exact %e)\n", approx_recourse, weighted_sum_of_recourse);
		}

		FREE_ARRAY_PTR(recourse_values);
	}

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

	/** allocate memory */
	SCIP_CALL(SCIPallocMemoryArray(scip, &vals, nvars_));

	/** get current solution */
	SCIP_CALL(SCIPgetSolVals(scip, sol, nvars_, vars_, vals));

	/** generate Benders cuts */
	generateCuts(nvars_, vals, where, &cs);

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
		DSPdebugMessage("naux_ %d nvars_ %d\n", naux_, nvars_);
		for (int j = nvars_ - naux_; j < nvars_; ++j)
		{
			if (cutrow.getMaxIndex() == j)
			{
				DSPdebugMessage("cutrow.getMaxIndex() = %d\n", j);
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
			SCIP_CALL(SCIPcreateEmptyRowConshdlr(scip, &row, conshdlr, "benders", rc->lb(), SCIPinfinity(scip),
					FALSE, /**< is row local? */
					FALSE, /**< is row modifiable? */
					BENDERS_CUT_REMOVABLE  /**< is row removable? can this be TRUE? */));

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
			DSPdebug(rc->print());

			/** flush all changes before adding cut */
			SCIP_CALL(SCIPflushRowExtensions(scip, row));

			DSPdebugMessage("efficacy %e isEfficatious %d\n", efficacy, isEfficacious);

			if (isEfficacious)
			{
				if (where == from_scip_enfolp)
				{
					/** add cut */
					SCIP_Bool infeasible;
					SCIP_CALL(SCIPaddRow(scip, row,
							BENDERS_CUT_FORCE, /**< force cut */
							&infeasible));

					if (infeasible || *result == SCIP_CUTOFF)
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
		else if (where != from_scip_enfolp && isEfficacious) {
			*result = SCIP_INFEASIBLE;
		}
	}
	DSPdebugMessage("sepaBenders: where %d result %d\n", where, *result);

	/** free memory */
	SCIPfreeMemoryArray(scip, &vals);

	return SCIP_OKAY;
}

void SCIPconshdlrBenders::setBdSub(BdSub * bdsub) {
	bdsub_ = bdsub;
	FREE_ARRAY_PTR(probability_);
	probability_ = new double [bdsub_->getNumSubprobs()];

	if (model_->isStochastic()) {
		// extract stochastic model
		TssModel* tss = dynamic_cast<TssModel*>(model_);
		CoinCopyN(tss->getProbability(), tss->getNumScenarios(), probability_);
	} else {
		CoinFillN(probability_, bdsub_->getNumSubprobs(), 1.0);
	}
}

/** set original variable pointers */
SCIP_RETCODE SCIPconshdlrBenders::setOriginalVariables(
		int nvars,        /**< number of original variables, including auxiliary variables */
		SCIP_Var ** vars, /**< original variables, including auxiliary variables */
		int         naux  /**< number of auxiliary variables */)
{
	nvars_ = nvars;
	SCIP_CALL(SCIPallocMemoryArray(scip_, &vars_, nvars_));
	for (int j = 0; j < nvars_; ++j)
		vars_[j] = vars[j];
	naux_ = naux;
	return SCIP_OKAY;
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

void SCIPconshdlrBenders::generateCuts(
		int size,      /**< [in] size of x */
		double * x,    /**< [in] master solution */
		int where,     /**< [in] where to be called */
		OsiCuts * cuts /**< [out] cuts generated */)
{
#define FREE_MEMORY                      \
	FREE_2D_ARRAY_PTR(nsubprobs, cutval) \
	FREE_ARRAY_PTR(cutrhs)

	assert(bdsub_);

	int nsubprobs = bdsub_->getNumSubprobs();
	double ** cutval = NULL;   /** dense cut coefficients for each subproblem */
	double *  cutrhs = NULL;   /** cut rhs for each subproblem */

	BGN_TRY_CATCH

	/** allocate memory */
	cutval = new double * [nsubprobs];
	cutrhs = new double [nsubprobs];

	/** generate cuts */
	// printf("x:\n");
	// DspMessage::printArray(nvars_, x);
	bdsub_->generateCuts(size, x, cutval, cutrhs);

	/** aggregate cuts */
	aggregateCuts(cutval, cutrhs, cuts);

#ifdef DSP_DEBUG
	printf("Generating cut at x:\n");
	DspMessage::printArray(size, x);
#endif
	DSPdebug(cuts->printCuts());

	END_TRY_CATCH(FREE_MEMORY)

	FREE_MEMORY

#undef FREE_MEMORY
}

void SCIPconshdlrBenders::aggregateCuts(
		double ** cutvec, /**< [in] cut vector */
		double *  cutrhs, /**< [in] cut right-hand side */
		OsiCuts * cuts    /**< [out] cuts generated */)
{
#define FREE_MEMORY                      \
	FREE_2D_ARRAY_PTR(naux_, aggval)     \
	FREE_ARRAY_PTR(aggrhs)

	int nsubprobs = bdsub_->getNumSubprobs();
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

	// TODO: The order of checking cuts should be the following, because of the feasibility cuts. 
	// The implementation can be vulnerable to any change to the cut generation algorithm.
	// This needs improved! -- Kibaek Kim 8/8/2018
	for (int i = nsubprobs - 1; i >= 0; --i)
	{
		DSPdebugMessage("bdsub_->getStatus(%d) = %d\n", i, bdsub_->getStatus(i));
		/** generate feasibility cut */
		if (bdsub_->getStatus(i) == DSP_STAT_PRIM_INFEASIBLE)
		{
#ifdef DSP_DEBUG
			printf("cutvec[%d]:\n", i);
			DspMessage::printArray(nvars_, cutvec[i]);
#endif
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
		if (bdsub_->getStatus(i) == DSP_STAT_NOT_SOLVED)
			break;

		if (bdsub_->getStatus(i) != DSP_STAT_OPTIMAL)
		{
			printf("Error: Subproblem %d returns unexpected status %d\n",
					bdsub_->getSubprobIndex(i), bdsub_->getStatus(i));
			for (int j = 0; j < cuts->sizeCuts(); ++j)
				delete cuts->rowCutPtr(i);
			cuts->dumpCuts();
			break;
		}

		int s = bdsub_->getSubprobIndex(i);         /**< subproblem index */
		int ind_aux = s % naux_;

		/** calculate weighted aggregation of cuts */
		for (int j = 0; j < nvars_; ++j)
			aggval[ind_aux][j] += cutvec[i][j] * probability_[s];
		aggrhs[ind_aux] += cutrhs[i] * probability_[s];
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

SCIP_RETCODE SCIPconshdlrBenders::evaluateRecourse(
		SCIP * scip,    /**< [in] scip pointer */
		SCIP_SOL * sol, /**< [in] solution to evaluate */
		double * values /**< [out] evaluated recourse values */) {
	/**< current solution */
	SCIP_Real * vals = NULL;

	/** allocate memory */
	SCIP_CALL(SCIPallocMemoryArray(scip, &vals, nvars_));

	/** get current solution */
	SCIP_CALL(SCIPgetSolVals(scip, sol, nvars_, vars_, vals));

	/** evaluate recourse */
	DSP_RTN_CODE ret = bdsub_->evaluateRecourse(vals, values);
	if (ret != DSP_RTN_OK) return SCIP_ERROR;

	/**
	 * TODO: the recourse values need to be de-scaled by the reference probability (for DRO).
	*/

	/** free memory */
	SCIPfreeMemoryArray(scip, &vals);

	return SCIP_OKAY;
}

void SCIPconshdlrBenders::computeProbability(
		const double* recourse /**< [in] recourse values */) {
	assert(model_->isStochastic());
	// TODO: find new probability distribution
}

SCIP_RETCODE SCIPconshdlrBenders::addNoGoodCut(
		SCIP * scip,             /**< [in] scip pointer */
		SCIP_CONSHDLR* conshdlr, /**< [in] constraint handler that creates the row */
		SCIP_RESULT * result     /**< [out] result */) {
	double rhs = 1.0;
	for (int j = 0; j < nvars_ - naux_; ++j) {
		rhs -= SCIPgetVarSol(scip, vars_[j]);
	}
	/** create empty row */
	SCIP_ROW * row = NULL;
	SCIP_CALL(SCIPcreateEmptyRowConshdlr(scip, &row, conshdlr, "No-good cut", rhs, SCIPinfinity(scip),
			FALSE, /**< is row local? */
			FALSE, /**< is row modifiable? */
			NOGOOD_CUT_REMOVABLE  /**< is row removable? can this be TRUE? */));

	/** cache the row extension and only flush them if the cut gets added */
	SCIP_CALL(SCIPcacheRowExtensions(scip, row));

	/** collect all non-zero coefficients */
	for (int j = 0; j < nvars_ - naux_; ++j) {
		if (SCIPisFeasZero(scip, SCIPgetVarSol(scip, vars_[j])))
			SCIP_CALL(SCIPaddVarToRow(scip, row, vars_[j], 1.0));
		else
			SCIP_CALL(SCIPaddVarToRow(scip, row, vars_[j], -1.0));
	}

	/** flush all changes before adding cut */
	SCIP_CALL(SCIPflushRowExtensions(scip, row));

	/** add cut */
	SCIP_Bool infeasible;
	SCIP_CALL(SCIPaddRow(scip, row, NOGOOD_CUT_FORCE, &infeasible));

	if (*result != SCIP_CUTOFF)
		*result = infeasible ? SCIP_CUTOFF : SCIP_SEPARATED;

	/** add cut to global pool */
	SCIP_CALL(SCIPaddPoolCut(scip, row));

	/** release the row */
	SCIP_CALL(SCIPreleaseRow(scip, &row));

	return SCIP_OKAY;
}

SCIP_RETCODE SCIPconshdlrBenders::addIntOptimalityCut(
		SCIP * scip,             /**< [in] scip pointer */
		SCIP_CONSHDLR* conshdlr, /**< [in] constraint handler that creates the row */
		double exact_recourse,   /**< [in] exact recourse value */
		double recourse_lb,  /**< [in] approximate recourse value */
		SCIP_RESULT * result     /**< [out] result */) {

	// compute right-hand side
	double ones = 0.0;
	for (int j = 0; j < nvars_ - naux_; ++j) {
		ones += SCIPgetVarSol(scip, vars_[j]);
	}
	double rhs = -ones*(exact_recourse-recourse_lb)+exact_recourse;

	/** create empty row */
	SCIP_ROW * row = NULL;
	SCIP_CALL(SCIPcreateEmptyRowConshdlr(scip, &row, conshdlr, "Integer Optimality cut", rhs, SCIPinfinity(scip),
			FALSE, /**< is row local? */
			FALSE, /**< is row modifiable? */
			INTOPT_CUT_REMOVABLE  /**< is row removable? can this be TRUE? */));

	/** cache the row extension and only flush them if the cut gets added */
	SCIP_CALL(SCIPcacheRowExtensions(scip, row));

	/** collect all non-zero coefficients */
	for (int j = 0; j < nvars_ - naux_; ++j) {
		if (SCIPisFeasZero(scip, SCIPgetVarSol(scip, vars_[j])))
			SCIP_CALL(SCIPaddVarToRow(scip, row, vars_[j], exact_recourse-recourse_lb));
		else
			SCIP_CALL(SCIPaddVarToRow(scip, row, vars_[j], recourse_lb-exact_recourse));
	}
	for (int j = 0; j < naux_; ++j) {
		SCIP_CALL(SCIPaddVarToRow(scip, row, vars_[nvars_-naux_+j], 1.0));
	}

	/** flush all changes before adding cut */
	SCIP_CALL(SCIPflushRowExtensions(scip, row));

	/** add cut */
	SCIP_Bool infeasible;
	SCIP_CALL(SCIPaddRow(scip, row, INTOPT_CUT_FORCE, &infeasible));

	if (*result != SCIP_CUTOFF)
		*result = infeasible ? SCIP_CUTOFF : SCIP_SEPARATED;

	/** add cut to global pool */
	SCIP_CALL(SCIPaddPoolCut(scip, row));

	/** release the row */
	SCIP_CALL(SCIPreleaseRow(scip, &row));

	return SCIP_OKAY;
}
