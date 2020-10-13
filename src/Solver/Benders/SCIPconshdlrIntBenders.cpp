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
							 true, /**< being in the initial LP? */
							 true,	/**< separated during LP process? */
							 true,	/**< enforced? */
							 true,	/**< checked for feasibility? */
							 true,	/**< propagate? */
							 false, /**< only locally valid? */
							 false, /**< modifiable? */
							 false, /**< is constraint subject to aging? */
							 false, /**< removable? */
							 false));

	return SCIP_OKAY;
}

/** default constructor */
SCIPconshdlrIntBenders::SCIPconshdlrIntBenders(SCIP *scip, const char *name, int sepapriority, int param_seps_solver)
	: SCIPconshdlrDrBenders(scip, name, sepapriority, param_seps_solver)
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
	SCIP_RESULT sepa_result = SCIP_FEASIBLE;
	SCIP_RESULT int_result = SCIP_FEASIBLE;
	SCIP_RESULT nogood_result = SCIP_FEASIBLE;
	
	/**
	 * This part computes the exact recourse function value at the current solution,
	 * 	adjusts the auxiliary variable values of the solution,
	 * 	and send it to SCIP as a primal solution.
	*/
	if (SCIPgetStage(scip) == SCIP_STAGE_SOLVING)
	{
		SCIP_Sol *sol;

		// get the current solution
		SCIP_CALL(SCIPcreateCurrentSol(scip, &sol, NULL));

		// Check whether this solution has been evaluated or not.
		if (check_binary_solution_pool(scip, sol) == false)
		{
			double *recourse_values = new double[model_->getNumSubproblems()];
		
			// evaluate the recourse value
			SCIP_CALL(evaluateRecourse(scip, sol, recourse_values));

			// compute weighted sum for DRO; otherwise, returns the current recourse
			if (model_->isDro())
				SCIP_CALL(computeProbability(recourse_values));

			// try solution after adjusting the auxiliary variable values
			SCIP_CALL(tryPrimalSolution(scip, sol, recourse_values));

			// try to add integer optimality cut
			SCIP_CALL(tryIntOptimalityCut(scip, conshdlr, sol, recourse_values, &int_result));
			DSPdebugMessage("----- scip_enfolp: integer Benders cut result %d\n", int_result);

			add_binary_solution_pool(scip, sol);
			
			FREE_ARRAY_PTR(recourse_values);
		}
		SCIP_CALL(SCIPfreeSol(scip, &sol));
	}

	SCIP_CALL(sepaBenders(scip, conshdlr, NULL, &sepa_result));
	DSPdebugMessage("----- scip_enfolp: Benders cut result %d stage %d\n", sepa_result, SCIPgetStage(scip));

	// Add no-good cut for pure binary
	addNoGoodCut(scip, conshdlr, &nogood_result);
	DSPdebugMessage("----- scip_enfolp: No-good cut result %d\n", nogood_result);
	
	// update the results
	*result = SCIP_FEASIBLE;
	if (int_result != SCIP_FEASIBLE)
		*result = int_result;
	if (sepa_result != SCIP_FEASIBLE && *result != SCIP_CUTOFF)
		*result = sepa_result;
	if (nogood_result != SCIP_FEASIBLE && *result != SCIP_CUTOFF)
		*result = nogood_result;

	return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for pseudo solutions */
SCIP_DECL_CONSENFOPS(SCIPconshdlrIntBenders::scip_enfops)
{
	*result = SCIP_FEASIBLE;
	SCIP_Sol *sol;

	// get the current solution
	SCIP_CALL(SCIPcreateCurrentSol(scip, &sol, NULL));

	/**
	 * This part computes the exact recourse function value at the current solution,
	 * 	compute the probability distribution,
	 * 	and check if the recourse approximate is valid.
	*/
	if (SCIPgetStage(scip) == SCIP_STAGE_SOLVING)
	{
		// Check whether this solution has been evaluated or not.
		if (check_binary_solution_pool(scip, sol) == false)
		{
			double *recourse_values = new double[model_->getNumSubproblems()];
		
			// evaluate the recourse value
			SCIP_CALL(evaluateRecourse(scip, sol, recourse_values));

			// compute weighted sum for DRO; otherwise, returns the current recourse
			if (model_->isDro())
				SCIP_CALL(computeProbability(recourse_values));

			double weighted_sum_of_recourse = compute_weighted_sum(recourse_values);
			double approx_recourse = compute_approximate_recourse(scip, sol, recourse_values);
			if (SCIPisLT(scip, approx_recourse, weighted_sum_of_recourse))
			{
				DSPdebugMessage("----- scip_enfops: rejects solution (approx %e, exact %e)\n", approx_recourse, weighted_sum_of_recourse);
				*result = SCIP_INFEASIBLE;
			}
			
			FREE_ARRAY_PTR(recourse_values);
		}
	}

	if (*result != SCIP_INFEASIBLE)
		SCIP_CALL(checkBenders(scip, conshdlr, sol, result));
	DSPdebugMessage("----- scip_enfops: result %d stage %d\n", *result, SCIPgetStage(scip));

	SCIP_CALL(SCIPfreeSol(scip, &sol));
	
	return SCIP_OKAY;
}

/** feasibility check method of constraint handler for primal solutions */
SCIP_DECL_CONSCHECK(SCIPconshdlrIntBenders::scip_check)
{
	*result = SCIP_FEASIBLE;

	/**
	 * This part computes the exact recourse function value at the current solution,
	 * 	compute the probability distribution,
	 * 	and check if the recourse approximate is valid.
	*/
	if (SCIPgetStage(scip) == SCIP_STAGE_SOLVING)
	{
		// Check whether this solution has been evaluated or not.
		if (check_binary_solution_pool(scip, sol) == false)
		{
			double *recourse_values = new double[model_->getNumSubproblems()];
		
			// evaluate the recourse value
			SCIP_CALL(evaluateRecourse(scip, sol, recourse_values));

			// compute weighted sum for DRO; otherwise, returns the current recourse
			if (model_->isDro())
				SCIP_CALL(computeProbability(recourse_values));

			double weighted_sum_of_recourse = compute_weighted_sum(recourse_values);
			double approx_recourse = compute_approximate_recourse(scip, sol, recourse_values);
			if (SCIPisLT(scip, approx_recourse, weighted_sum_of_recourse))
			{
				DSPdebugMessage("----- scip_check: rejects solution (approx %e, exact %e)\n", approx_recourse, weighted_sum_of_recourse);
				*result = SCIP_INFEASIBLE;
			}
			// else
			// {
			// 	printf("----- scip_check: accepts solution (approx %e, exact %e)\n", approx_recourse, weighted_sum_of_recourse);
			// 	SCIPprintSol(scip, sol, NULL, FALSE);
			// }

			FREE_ARRAY_PTR(recourse_values);
		}
		else
			*result = SCIP_INFEASIBLE;
	}

	if (*result != SCIP_INFEASIBLE)
		SCIP_CALL(checkBenders(scip, conshdlr, sol, result));
	// DSPdebugMessage("----- scip_check: result %d stage %d\n", *result, SCIPgetStage(scip));

	return SCIP_OKAY;
}

SCIP_DECL_CONSSEPALP(SCIPconshdlrIntBenders::scip_sepalp)
{
	*result = SCIP_DIDNOTFIND;
	if (model_->isDro())
		SCIP_CALL(sepaDrBenders(scip, conshdlr, NULL, result));
	else
		SCIP_CALL(sepaBenders(scip, conshdlr, NULL, result));
	if (*result == SCIP_DIDNOTRUN)
		*result = SCIP_DIDNOTFIND;
	DSPdebugMessage("scip_sepalp: results %d stage %d\n", *result, SCIPgetStage(scip));
	return SCIP_OKAY;
}

SCIP_DECL_CONSSEPASOL(SCIPconshdlrIntBenders::scip_sepasol)
{
	*result = SCIP_DIDNOTFIND;
	if (model_->isDro())
		SCIP_CALL(sepaDrBenders(scip, conshdlr, NULL, result));
	else
		SCIP_CALL(sepaBenders(scip, conshdlr, NULL, result));
	if (*result == SCIP_DIDNOTRUN)
		*result = SCIP_DIDNOTFIND;
	DSPdebugMessage("scip_sepasol results in %d stage %d\n", *result, SCIPgetStage(scip));
	return SCIP_OKAY;
}

double SCIPconshdlrIntBenders::compute_weighted_sum(
	const double *recourse_values /**< [in] recourse values */)
{
	double v = 0.;
	for (int j = 0; j < model_->getNumSubproblems(); ++j)
	{
		// DSPdebugMessage("----- scip_enfops: exact recourse value [%d] %e with probability %e\n", j, recourse_values[j], probability_[j]);
		v += recourse_values[j] * probability_[j];
	}
	return v;
}

double SCIPconshdlrIntBenders::compute_approximate_recourse(
	SCIP *scip,	   /**< [in] scip pointer */
	SCIP_SOL *sol, /**< [in] solution to evaluate */
	const double *recourse_values /**< [in] recourse values */)
{
	double v = 0.0;
	for (int j = 0; j < naux_; ++j)
		v += SCIPgetSolVal(scip, sol, vars_[nvars_ - naux_ + j]);
	return v;
}

SCIP_RETCODE SCIPconshdlrIntBenders::tryPrimalSolution(
		SCIP *scip,	   /**< [in] scip pointer */
		SCIP_SOL *sol, /**< [in] solution to evaluate */
		const double *recourse_values /**< [in] recourse values */)
{
	SCIP_Bool stored;
	
	// set a feasible value of the auxiliary variable
	for (int j = 0; j < naux_; ++j)
	{
		double weighted_recourse = 0.0;
		for (int k = j; k < model_->getNumSubproblems(); k += naux_)
		{
			weighted_recourse += recourse_values[k] * probability_[k];
			// DSPdebugMessage("---- s = %d, recourse_values = %e, probability_ %e\n", k, recourse_values[k], probability_[k]);
		}
		SCIP_CALL(SCIPsetSolVal(scip, sol, vars_[nvars_ - naux_ + j], weighted_recourse));
	}
	// SCIPprintSol(scip, sol, NULL, FALSE);

	// set the solution as a primal feasible solution
	double sense = SCIPgetObjsense(scip) == SCIP_OBJSENSE_MINIMIZE ? 1.0 : -1.0;
	if (SCIPisLT(scip, (SCIPsolGetOrigObj(sol) - SCIPgetPrimalbound(scip)) * sense, 0.0))
	{
		SCIP_CALL(SCIPtrySol(scip, sol, FALSE, FALSE, FALSE, FALSE, FALSE, &stored));
		DSPdebugMessage("----- is solution (obj: %e) stored? %s\n", SCIPsolGetOrigObj(sol), stored ? "yes" : "no");
	}
	
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

SCIP_RETCODE SCIPconshdlrIntBenders::tryIntOptimalityCut(
	SCIP *scip,				 /**< [in] scip pointer */
	SCIP_CONSHDLR *conshdlr, /**< [in] constraint handler that creates the row */
	SCIP_SOL *sol, /**< [in] solution to evaluate */
	const double *recourse_values, /**< [in] recourse values */
	SCIP_RESULT *result /**< [out] result */)
{
	double weighted_sum_of_recourse = compute_weighted_sum(recourse_values);
	double approx_recourse = compute_approximate_recourse(scip, sol, recourse_values);
	if (SCIPisLT(scip, approx_recourse, weighted_sum_of_recourse))
	{
		addIntOptimalityCut(scip, conshdlr, weighted_sum_of_recourse, result);
		DSPdebugMessage("---- integer optimality cut result %d\n", *result);
	}
	return SCIP_OKAY;
}

SCIP_RETCODE SCIPconshdlrIntBenders::addIntOptimalityCut(
	SCIP *scip,				 /**< [in] scip pointer */
	SCIP_CONSHDLR *conshdlr, /**< [in] constraint handler that creates the row */
	double exact_recourse,	 /**< [in] exact recourse value */
	SCIP_RESULT *result /**< [out] result */)
{

	double stime = CoinGetTimeOfDay();

	// compute the lower bound of the recourse function
	double recourse_lb = SCIPgetDualbound(scip);
	for (int j = 0; j < nvars_ - naux_; ++j)
		recourse_lb -= CoinMax(0.0, SCIPvarGetObj(vars_[j]));

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
	SCIP_SOL *sol /**< [in] solution to evaluate */)
{
	bool exist = false;
	// return exist;
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
	return exist;
}

void SCIPconshdlrIntBenders::add_binary_solution_pool(
	SCIP *scip,	   /**< [in] scip pointer */
	SCIP_SOL *sol /**< [in] solution to evaluate */)
{
	vector<int> solvec;
	solvec.reserve(nvars_ - naux_);
	for (int j = 0; j < nvars_ - naux_; ++j)
		if (SCIPisZero(scip, SCIPgetSolVal(scip, sol, vars_[j]) - 1.0))
			solvec.push_back(j);
	binary_solution_pool_.push_back(solvec);
}