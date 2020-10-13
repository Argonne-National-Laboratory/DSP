/*
 * SCIPconshdlrDrBenders.cpp
 */

// #define DSP_DEBUG
// #define BENDERS_PROFILE

#include "OsiCuts.hpp"

/** DSP */
#include "Utility/DspMessage.h"
#include "Model/TssModel.h"
#include "Solver/Benders/SCIPconshdlrDrBenders.h"
#include "SolverInterface/DspOsiScip.h"
#include "SolverInterface/DspOsiClp.h"
#include "SolverInterface/DspOsiCpx.h"
#include "SolverInterface/DspOsiGrb.h"

#define BENDERS_CUT_REMOVABLE FALSE
#define BENDERS_CUT_FORCE TRUE

/** constraint data for Benders cuts */
struct SCIP_ConsData
{
	/** nothing */
};

/** creates and captures a distributionally robust Benders constraint */
SCIP_RETCODE SCIPcreateConsDrBenders(
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
							 true,	/**< being in the initial LP? */
							 true,	/**< separated during LP process? */
							 true,	/**< enforced? */
							 true,	/**< checked for feasibility? */
							 true,	/**< propagate? */
							 false, /**< only locally valid? */
							 false, /**< modifiable? */
							 false, /**< is constraint subject to aging? */
							 false,	/**< removable? */
							 false));
	return SCIP_OKAY;
}

/** default constructor */
SCIPconshdlrDrBenders::SCIPconshdlrDrBenders(SCIP *scip, const char *name, int sepapriority, int param_seps_solver)
	: SCIPconshdlrBenders(scip, name, sepapriority)
{
	switch (param_seps_solver)
	{
	case OsiCpx:
#ifdef DSP_HAS_CPX
		drosi_ = new DspOsiCpx();
#else
		throw CoinError("Cplex is not available.", "SCIPconshdlrDrBenders", "SCIPconshdlrDrBenders");
#endif
		break;
	case OsiGrb:
#ifdef DSP_HAS_GRB
		drosi_ = new DspOsiGrb();
#else
		throw CoinError("Gurobi is not available.", "SCIPconshdlrDrBenders", "SCIPconshdlrDrBenders");
#endif
		break;
	case OsiClp:
		drosi_ = new DspOsiClp();
		break;
	default:
		char coinmsg[128];
		sprintf(coinmsg, "Invalid parameter value (solver = %d)", param_seps_solver);
		throw CoinError(coinmsg, "SCIPconshdlrDrBenders", "SCIPconshdlrDrBenders");
		break;
	}
#ifdef BENDERS_PROFILE
	/** initialize statistics */
	names_statistics_.push_back("sepaBenders");
	names_statistics_.push_back("checkBenders");
	names_statistics_.push_back("evaluateRecourse");
	for (unsigned i = 0; i < names_statistics_.size(); ++i)
	{
		time_statistics_[names_statistics_[i]] = 0.0;
		count_statistics_[names_statistics_[i]] = 0;
	}
#endif
}

SCIPconshdlrDrBenders::~SCIPconshdlrDrBenders()
{
#ifdef BENDERS_PROFILE
	write_statistics();
#endif
	FREE_PTR(drosi_);
}

/** constraint enforcing method of constraint handler for LP solutions */
SCIP_DECL_CONSENFOLP(SCIPconshdlrDrBenders::scip_enfolp)
{
	*result = SCIP_FEASIBLE;
	DSPdebugMessage("scip_enfolp\n");
	SCIP_CALL(sepaDrBenders(scip, conshdlr, NULL, result));
	DSPdebugMessage("scip_enfolp: result %d, stage %d\n", *result, SCIPgetStage(scip));

	return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for pseudo solutions */
SCIP_DECL_CONSENFOPS(SCIPconshdlrDrBenders::scip_enfops)
{
	*result = SCIP_FEASIBLE;
	SCIP_CALL(checkDrBenders(scip, conshdlr, NULL, result));
	if (*result == SCIP_SEPARATED)
		*result = SCIP_INFEASIBLE;
	return SCIP_OKAY;
}

/** feasibility check method of constraint handler for primal solutions */
SCIP_DECL_CONSCHECK(SCIPconshdlrDrBenders::scip_check)
{
	*result = SCIP_FEASIBLE;
	DSPdebugMessage("scip_check\n");
	SCIP_CALL(checkDrBenders(scip, conshdlr, sol, result));
	DSPdebugMessage("scip_check results in %d stage(%d)\n", *result, SCIPgetStage(scip));
	return SCIP_OKAY;
}

SCIP_DECL_CONSSEPALP(SCIPconshdlrDrBenders::scip_sepalp)
{
	*result = SCIP_DIDNOTFIND;
	DSPdebugMessage("scip_sepalp\n");
	SCIP_CALL(sepaDrBenders(scip, conshdlr, NULL, result));
	DSPdebugMessage("scip_sepalp results in %d stage %d\n", *result, SCIPgetStage(scip));
	return SCIP_OKAY;
}

SCIP_DECL_CONSSEPASOL(SCIPconshdlrDrBenders::scip_sepasol)
{
	*result = SCIP_DIDNOTFIND;
	SCIP_CALL(sepaDrBenders(scip, conshdlr, NULL, result));
	DSPdebugMessage("scip_sepasol results in %d stage %d\n", *result, SCIPgetStage(scip));
	return SCIP_OKAY;
}

void SCIPconshdlrDrBenders::setDecModel(DecModel *model)
{
	model_ = model;
	create_distsepa_problem();
}

SCIP_RETCODE SCIPconshdlrDrBenders::evaluateRecourse(
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

void SCIPconshdlrDrBenders::create_distsepa_problem()
{
#define FREE_MEMORY      \
	FREE_PTR(mat)        \
	FREE_ARRAY_PTR(clbd) \
	FREE_ARRAY_PTR(cubd) \
	FREE_ARRAY_PTR(obj)  \
	FREE_ARRAY_PTR(rlbd) \
	FREE_ARRAY_PTR(rubd)

	CoinPackedMatrix *mat = NULL;
	double *clbd = NULL;
	double *cubd = NULL;
	double *obj = NULL;
	double *rlbd = NULL;
	double *rubd = NULL;

	BGN_TRY_CATCH

	int nscen = model_->getNumSubproblems();
	int nrefs = model_->getNumReferences();
	int ncols = nscen + nscen * nrefs;
	int nrows = 1 + nscen + nrefs + 1;

	clbd = new double[ncols];
	cubd = new double[ncols];
	obj = new double[ncols];
	rlbd = new double[nrows];
	rubd = new double[nrows];

	CoinZeroN(clbd, ncols);
	CoinFillN(cubd, ncols, COIN_DBL_MAX);
	CoinZeroN(obj, ncols);

	// rhss
	rlbd[0] = -COIN_DBL_MAX;
	rubd[0] = model_->getWassersteinSize();
	CoinZeroN(rlbd + 1, nscen);
	CoinZeroN(rubd + 1, nscen);
	for (int j = 0; j < nrefs; ++j)
	{
		rlbd[1 + nscen + j] = model_->getReferenceProbability(j);
		rubd[1 + nscen + j] = model_->getReferenceProbability(j);
	}
	rlbd[nrows - 1] = 1.0;
	rubd[nrows - 1] = 1.0;

	mat = new CoinPackedMatrix(false, 0, 0);
	mat->setDimensions(0, ncols);

	vector<int> vind;
	vector<double> velem;
	vind.reserve(nscen * nrefs);
	velem.reserve(nscen * nrefs);

	// Wasserstein ball constraint
	for (int i = 0; i < nscen; ++i)
		for (int j = 0; j < nrefs; ++j)
		{
			vind.push_back(nscen + nrefs * i + j);
			velem.push_back(model_->getWassersteinDist(j, i));
		}
	mat->appendRow(nrefs * nscen, vind.data(), velem.data());

	// marginal probability of scenario
	velem.clear();
	velem.push_back(-1.0);
	for (int j = 0; j < nrefs; ++j)
		velem.push_back(1.0);
	for (int i = 0; i < nscen; ++i)
	{
		vind.clear();
		vind.push_back(i);
		for (int j = 0; j < nrefs; ++j)
			vind.push_back(nscen + nrefs * i + j);
		mat->appendRow(nrefs + 1, vind.data(), velem.data());
	}

	// marginal probability of reference
	velem.clear();
	for (int i = 0; i < nscen; ++i)
		velem.push_back(1.0);
	for (int j = 0; j < nrefs; ++j)
	{
		vind.clear();
		for (int i = 0; i < nscen; ++i)
			vind.push_back(nscen + nrefs * i + j);
		mat->appendRow(nscen, vind.data(), velem.data());
	}

	// probability summing up to one
	vind.clear();
	velem.clear();
	for (int i = 0; i < nscen; ++i)
	{
		vind.push_back(i);
		velem.push_back(1.0);
	}
	mat->appendRow(nscen, vind.data(), velem.data());
	DSPdebug(mat->verifyMtx(4));

	drosi_->si_->loadProblem(*mat, clbd, cubd, obj, rlbd, rubd);
	drosi_->si_->setObjSense(-1.0);
	drosi_->setLogLevel(0);

	END_TRY_CATCH(FREE_MEMORY)

	FREE_MEMORY;

#undef FREE_MEMORY
}

SCIP_RETCODE SCIPconshdlrDrBenders::computeProbability(
	const double *recourse /**< [in] recourse values */)
{
	assert(isStochastic());
	assert(drosi_);
	assert(drosi_->si_);
	assert(drosi_->si_->getNumCols() >= model_->getNumSubproblems());
	// change objective function
	for (int j = 0; j < model_->getNumSubproblems(); ++j)
		drosi_->si_->setObjCoeff(j, recourse[j]);
	DSPdebugMessage("updated distribution separation objective.\n");

	// solve the distribution separation problem
	drosi_->solve();
	DSPdebugMessage("distribution separation status: %d\n", drosi_->status());

	if (drosi_->status() == DSP_STAT_OPTIMAL)
	{
		for (int j = 0; j < model_->getNumSubproblems(); ++j)
		{
			probability_[j] = drosi_->si_->getColSolution()[j];
			DSPdebugMessage("probability_[%d] = %e\n", j, probability_[j]);
		}
	}
	else
	{
		printf("Unexpected status: %d\n", drosi_->status());
		drosi_->si_->writeMps("drsepa");
		return SCIP_ERROR;
	}

	return SCIP_OKAY;
}

SCIP_RETCODE SCIPconshdlrDrBenders::sepaDrBenders(
	SCIP *scip, /**< [in] scip pointer */
	SCIP_CONSHDLR *conshdlr,
	SCIP_SOL *sol, /**< [in] solution to evaluate */
	SCIP_RESULT *result /**< [out] result */)
{
	/**
	 * This part computes the exact recourse function value at the current solution
	 * 	and adjusts the probability distribution.
	*/
	if (SCIPgetStage(scip) == SCIP_STAGE_SOLVING)
	{
		SCIP_Sol *tmpsol;

		// get the current solution
		if (sol == NULL)
			SCIP_CALL(SCIPcreateCurrentSol(scip, &tmpsol, NULL));
		else
			tmpsol = sol;

		double *recourse_values = new double[model_->getNumSubproblems()];

		// evaluate the recourse value
		SCIP_CALL(evaluateRecourse(scip, tmpsol, recourse_values));

		// compute weighted sum for DRO; otherwise, returns the current recourse
		SCIP_CALL(computeProbability(recourse_values));

		// free memeory
		FREE_ARRAY_PTR(recourse_values);

		// free solution memeory
		if (sol == NULL)
			SCIP_CALL(SCIPfreeSol(scip, &tmpsol));
		else
			tmpsol = NULL;
	}

	SCIP_CALL(sepaBenders(scip, conshdlr, sol, result));

	return SCIP_OKAY;
}

SCIP_RETCODE SCIPconshdlrDrBenders::checkDrBenders(
	SCIP *scip, /**< [in] scip pointer */
	SCIP_CONSHDLR *conshdlr,
	SCIP_SOL *sol, /**< [in] solution to evaluate */
	SCIP_RESULT *result /**< [out] result */)
{
	/**
	 * This part computes the exact recourse function value at the current solution
	 * 	and adjusts the probability distribution.
	*/
	if (SCIPgetStage(scip) == SCIP_STAGE_SOLVING)
	{
		SCIP_Sol *tmpsol;

		// get the current solution
		if (sol == NULL)
			SCIP_CALL(SCIPcreateCurrentSol(scip, &tmpsol, NULL));
		else
			tmpsol = sol;

		double *recourse_values = new double[model_->getNumSubproblems()];

		// evaluate the recourse value
		SCIP_CALL(evaluateRecourse(scip, tmpsol, recourse_values));

		// compute weighted sum for DRO; otherwise, returns the current recourse
		SCIP_CALL(computeProbability(recourse_values));

		// free memeory
		FREE_ARRAY_PTR(recourse_values);

		// free solution memeory
		if (sol == NULL)
			SCIP_CALL(SCIPfreeSol(scip, &tmpsol));
		else
			tmpsol = NULL;
	}

	SCIP_CALL(checkBenders(scip, conshdlr, sol, result));

	return SCIP_OKAY;
}
