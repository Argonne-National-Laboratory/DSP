/*
 * SCIPconshdlrIntBenders.cpp
 */

// #define DSP_DEBUG
// #define BENDERS_PROFILE

/** DSP */
#include "Utility/DspMessage.h"
#include "Model/TssModel.h"
#include "Solver/Benders/SCIPconshdlrIntBenders.h"

#define INTOPT_CUT_REMOVABLE FALSE
#define INTOPT_CUT_FORCE TRUE
#define NOGOOD_CUT_REMOVABLE FALSE
#define NOGOOD_CUT_FORCE TRUE

/** constraint data for Benders cuts */
struct SCIP_ConsData
{
	/** nothing */
};

/** creates and captures a Benders constraint */
SCIP_RETCODE SCIPcreateConsIntBenders(
	SCIP *scip,
	SCIP_CONS **cons,
	const char *name)
{
	SCIP_CONSHDLR *conshdlr = NULL;
	SCIP_CONSDATA *consdata = NULL;

	/* find the Benders constraint handler */
	conshdlr = SCIPfindConshdlr(scip, name);
	assert(conshdlr);

	/* create constraint */
	SCIP_CALL(SCIPcreateCons(scip, cons, name, conshdlr, consdata,
							 false, /**< being in the initial LP? */
							 true,	/**< separated during LP process? */
							 true,	/**< enforced? */
							 true,	/**< checked for feasibility? */
							 true,	/**< propagate? */
							 false, /**< only locally valid? */
							 false, /**< modifiable? */
							 false, /**< is constraint subject to aging? */
							 true,	/**< removable? */
							 false));

	return SCIP_OKAY;
}

/** default constructor */
SCIPconshdlrIntBenders::SCIPconshdlrIntBenders(
	SCIP *scip,
	const char *name,
	int sepapriority)
	: SCIPconshdlrBenders(scip, name, sepapriority)
{
#ifdef BENDERS_PROFILE
	/** initialize statistics */
	names_statistics_.push_back("sepaBenders");
	names_statistics_.push_back("checkBenders");
	names_statistics_.push_back("addIntOptimalityCut");
	names_statistics_.push_back("addNoGoodCut");
	names_statistics_.push_back("evaluateRecourse");
	for (unsigned i = 0; i < names_statistics_.size(); ++i)
	{
		time_statistics_[names_statistics_[i]] = 0.0;
		count_statistics_[names_statistics_[i]] = 0;
	}
#endif
}

SCIPconshdlrIntBenders::~SCIPconshdlrIntBenders()
{
#ifdef BENDERS_PROFILE
	write_statistics();
#endif
}

/** constraint enforcing method of constraint handler for LP solutions */
SCIP_DECL_CONSENFOLP(SCIPconshdlrIntBenders::scip_enfolp)
{
#define FREE_MEMORY \
	FREE_ARRAY_PTR(recourse_values)

	double *recourse_values = NULL;
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
	if (SCIPgetStage(scip) == SCIP_STAGE_SOLVING)
	{
		SCIP_Bool stored;
		SCIP_Sol *sol;
		nsubs = model_->getNumSubproblems();

		recourse_values = new double[nsubs];

		// get the current solution
		SCIP_CALL(SCIPcreateCurrentSol(scip, &sol, NULL));

		// Check whether this solution has been evaluated or not.
		if (check_binary_solution_pool(scip, sol, true) == false)
		{
			// evaluate the recourse value
			SCIP_CALL(evaluateRecourse(scip, sol, recourse_values));

			// compute weighted sum for DRO; otherwise, returns the current recourse
			computeProbability(recourse_values);
			for (int j = 0; j < nsubs; ++j)
			{
				DSPdebugMessage("----- scip_enfolp: exact recourse value [%d] %e with probability %e\n", j, recourse_values[j], probability_[j]);
				weighted_sum_of_recourse += recourse_values[j] * probability_[j];
			}

			// set a feasible value of the auxiliary variable
			for (int j = 0; j < naux_; ++j)
			{
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
			// SCIPprintSol(scip, sol, NULL, FALSE);

			// set the solution as a primal feasible solution
			double sense = SCIPgetObjsense(scip) == SCIP_OBJSENSE_MINIMIZE ? 1.0 : -1.0;
			if (SCIPisLT(scip, (SCIPsolGetOrigObj(sol) - SCIPgetPrimalbound(scip)) * sense, 0.0))
			{
				SCIP_CALL(SCIPtrySol(scip, sol, FALSE, FALSE, FALSE, FALSE, FALSE, &stored));
				DSPdebugMessage("----- scip_enfolp: is solution (obj: %e) stored? %s\n", SCIPsolGetOrigObj(sol), stored ? "yes" : "no");
			}
		}
		SCIP_CALL(SCIPfreeSol(scip, &sol));
	}

	*result = SCIP_FEASIBLE;
	/**
	 * TODO: DRO needs to pass probability_ to sepaBenders
	 */
	SCIP_CALL(sepaBenders(scip, conshdlr, NULL, result));
	DSPdebugMessage("scip_enfolp results in %d stage %d\n", *result, SCIPgetStage(scip));

	if (SCIPgetStage(scip) == SCIP_STAGE_SOLVING &&
		approx_recourse < weighted_sum_of_recourse)
	{
		/** integer Benders cut */
		double recourse_lb = SCIPgetDualbound(scip);
		for (int j = 0; j < nvars_ - naux_; ++j)
		{
			recourse_lb -= CoinMax(0.0, SCIPvarGetObj(vars_[j]));
		}
		addIntOptimalityCut(scip, conshdlr, weighted_sum_of_recourse, recourse_lb, result);
		DSPdebugMessage("----- scip_enfolp: integer optimality cut result %d\n", *result);
	}

	// Add no-good cut for pure binary
	addNoGoodCut(scip, conshdlr, result);
	DSPdebugMessage("----- scip_enfolp: No-good cut result %d\n", *result);

	FREE_MEMORY;

	return SCIP_OKAY;
#undef FREE_MEMORY
}

/** constraint enforcing method of constraint handler for pseudo solutions */
SCIP_DECL_CONSENFOPS(SCIPconshdlrIntBenders::scip_enfops)
{
	*result = SCIP_FEASIBLE;
	SCIP_CALL(sepaBenders(scip, conshdlr, NULL, result));
	if (*result == SCIP_SEPARATED)
		*result = SCIP_INFEASIBLE;
	return SCIP_OKAY;
}

/** feasibility check method of constraint handler for primal solutions */
SCIP_DECL_CONSCHECK(SCIPconshdlrIntBenders::scip_check)
{
	*result = SCIP_FEASIBLE;

	SCIP_CALL(checkBenders(scip, conshdlr, sol, result));
	DSPdebugMessage("scip_check results in %d stage(%d)\n", *result, SCIPgetStage(scip));

	if (SCIPgetStage(scip) == SCIP_STAGE_SOLVING &&
		*result == SCIP_FEASIBLE &&
		check_binary_solution_pool(scip, sol) == false)
	{
		int nsubs = model_->getNumSubproblems();
		double weighted_sum_of_recourse = 0.0;
		double approx_recourse = 0.0;

		double *recourse_values = new double[nsubs];

		// evaluate the recourse value
		SCIP_CALL(evaluateRecourse(scip, sol, recourse_values));

		// compute weighted sum for DRO; otherwise, returns the current recourse
		computeProbability(recourse_values);
		for (int j = 0; j < nsubs; ++j)
		{
			DSPdebugMessage("----- scip_check: exact recourse value [%d] %e with probability %e\n", j, recourse_values[j], probability_[j]);
			weighted_sum_of_recourse += recourse_values[j] * probability_[j];
		}

		// set a feasible value of the auxiliary variable
		for (int j = 0; j < naux_; ++j)
		{
			approx_recourse += SCIPgetSolVal(scip, sol, vars_[nvars_ - naux_ + j]);
		}
		// DSPdebugMessage("----- scip_check:\n");
		// SCIPprintSol(scip, sol, NULL, FALSE);
		// DSPdebugMessage("----- scip_check: recourse values (approx %e, exact %e)\n", approx_recourse, weighted_sum_of_recourse);

		if (approx_recourse < weighted_sum_of_recourse)
		{
			*result = SCIP_INFEASIBLE;
			DSPdebugMessage("----- scip_check: rejects solution (approx %e, exact %e)\n", approx_recourse, weighted_sum_of_recourse);
		}

		FREE_ARRAY_PTR(recourse_values);
	}

	return SCIP_OKAY;
}

SCIP_RETCODE SCIPconshdlrIntBenders::evaluateRecourse(
	SCIP *scip,	   /**< [in] scip pointer */
	SCIP_SOL *sol, /**< [in] solution to evaluate */
	double *values /**< [out] evaluated recourse values */)
{

	double stime = CoinGetTimeOfDay();

	/**< current solution */
	SCIP_Real *vals = NULL;

	/** allocate memory */
	SCIP_CALL(SCIPallocMemoryArray(scip, &vals, nvars_));

	/** get current solution */
	SCIP_CALL(SCIPgetSolVals(scip, sol, nvars_, vars_, vals));

	/** evaluate recourse */
	DSP_RTN_CODE ret = bdsub_->evaluateRecourse(vals, values);
	if (ret != DSP_RTN_OK)
		return SCIP_ERROR;

	/** free memory */
	SCIPfreeMemoryArray(scip, &vals);

#ifdef BENDERS_PROFILE
	time_statistics_["evaluateRecourse"] += CoinGetTimeOfDay() - stime;
	count_statistics_["evaluateRecourse"]++;
#endif

	return SCIP_OKAY;
}

SCIP_RETCODE SCIPconshdlrIntBenders::addNoGoodCut(
	SCIP *scip,				 /**< [in] scip pointer */
	SCIP_CONSHDLR *conshdlr, /**< [in] constraint handler that creates the row */
	SCIP_RESULT *result /**< [out] result */)
{

	double stime = CoinGetTimeOfDay();

	double rhs = 1.0;
	for (int j = 0; j < nvars_ - naux_; ++j)
	{
		rhs -= SCIPgetVarSol(scip, vars_[j]);
	}
	/** create empty row */
	SCIP_ROW *row = NULL;
	SCIP_CALL(SCIPcreateEmptyRowConshdlr(scip, &row, conshdlr, "No-good cut", rhs, SCIPinfinity(scip),
										 FALSE, /**< is row local? */
										 FALSE, /**< is row modifiable? */
										 NOGOOD_CUT_REMOVABLE /**< is row removable? can this be TRUE? */));

	/** cache the row extension and only flush them if the cut gets added */
	SCIP_CALL(SCIPcacheRowExtensions(scip, row));

	/** collect all non-zero coefficients */
	for (int j = 0; j < nvars_ - naux_; ++j)
	{
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

#ifdef BENDERS_PROFILE
	time_statistics_["addNoGoodCut"] += CoinGetTimeOfDay() - stime;
	count_statistics_["addNoGoodCut"]++;
#endif

	return SCIP_OKAY;
}

SCIP_RETCODE SCIPconshdlrIntBenders::addIntOptimalityCut(
	SCIP *scip,				 /**< [in] scip pointer */
	SCIP_CONSHDLR *conshdlr, /**< [in] constraint handler that creates the row */
	double exact_recourse,	 /**< [in] exact recourse value */
	double recourse_lb,		 /**< [in] approximate recourse value */
	SCIP_RESULT *result /**< [out] result */)
{

	double stime = CoinGetTimeOfDay();

	// compute right-hand side
	double ones = 0.0;
	for (int j = 0; j < nvars_ - naux_; ++j)
	{
		ones += SCIPgetVarSol(scip, vars_[j]);
	}
	double rhs = -ones * (exact_recourse - recourse_lb) + exact_recourse;

	/** create empty row */
	SCIP_ROW *row = NULL;
	SCIP_CALL(SCIPcreateEmptyRowConshdlr(scip, &row, conshdlr, "Integer Optimality cut", rhs, SCIPinfinity(scip),
										 FALSE, /**< is row local? */
										 FALSE, /**< is row modifiable? */
										 INTOPT_CUT_REMOVABLE /**< is row removable? can this be TRUE? */));

	/** cache the row extension and only flush them if the cut gets added */
	SCIP_CALL(SCIPcacheRowExtensions(scip, row));

	/** collect all non-zero coefficients */
	for (int j = 0; j < nvars_ - naux_; ++j)
	{
		if (SCIPisFeasZero(scip, SCIPgetVarSol(scip, vars_[j])))
			SCIP_CALL(SCIPaddVarToRow(scip, row, vars_[j], exact_recourse - recourse_lb));
		else
			SCIP_CALL(SCIPaddVarToRow(scip, row, vars_[j], recourse_lb - exact_recourse));
	}
	for (int j = 0; j < naux_; ++j)
	{
		SCIP_CALL(SCIPaddVarToRow(scip, row, vars_[nvars_ - naux_ + j], 1.0));
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

#ifdef BENDERS_PROFILE
	time_statistics_["addIntOptimalityCut"] += CoinGetTimeOfDay() - stime;
	count_statistics_["addIntOptimalityCut"]++;
#endif

	return SCIP_OKAY;
}

bool SCIPconshdlrIntBenders::check_binary_solution_pool(
	SCIP *scip,	   /**< [in] scip pointer */
	SCIP_SOL *sol, /**< [in] solution to evaluate */
	bool append /**< [in] whether sol is appended to the pool */)
{
	bool exist = false;
	vector<int> solvec;
	solvec.reserve(nvars_ - naux_);
	for (int j = 0; j < nvars_ - naux_; ++j)
		if (SCIPisZero(scip, SCIPgetSolVal(scip, sol, vars_[j]) - 1.0))
			solvec.push_back(j);
	for (unsigned i = 0; i < binary_solution_pool_.size(); ++i)
		if (solvec == binary_solution_pool_[i])
		{
			exist = true;
			break;
		}
	if (!exist && append)
		binary_solution_pool_.push_back(solvec);
	return exist;
}