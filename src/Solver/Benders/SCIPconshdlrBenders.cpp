/*
 * SCIPconshdlrBenders.cpp
 *
 *  Created on: Dec 8, 2014
 *      Author: kibaekkim
 */

// #define DSP_DEBUG
// #define BENDERS_PROFILE

#include "OsiCuts.hpp"

/** DSP */
#include "Utility/DspMessage.h"
#include "Model/TssModel.h"
#include "Solver/Benders/SCIPconshdlrBenders.h"

#define BENDERS_CUT_REMOVABLE FALSE
#define BENDERS_CUT_FORCE TRUE

/** constraint data for Benders cuts */
struct SCIP_ConsData
{
	/** nothing */
};

/** creates and captures a Benders constraint */
SCIP_RETCODE SCIPcreateConsBenders(
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
							 false,	/**< removable? */
							 false));
	return SCIP_OKAY;
}

/** default constructor */
SCIPconshdlrBenders::SCIPconshdlrBenders(
	SCIP *scip,
	const char *name,
	int sepapriority)
	: ObjConshdlr(scip, name, "Benders cuts",
				  1, /**< priority of the constraint handler for separation */
				  sepapriority, /**< priority of the constraint handler for constraint enforcing */
				  sepapriority, /**< priority of the constraint handler for checking infeasibility (and propagation) */
				  0,			/**< run separator only at root */
				  1,			/**< frequency for propagating domains; zero means only preprocessing propagation */
				  1,			/**< always use all constraints (no aging) */
				  0,			/**< disable the preprocessing callback of the constraint handler */
				  TRUE,			/**< delay separation method */
				  FALSE,		/**< do not delay propatation method */
				  TRUE,			/**< skip constraint handler even if no constraints are available. */
#if SCIP_VERSION < 320
				  SCIP_PROPTIMING_BEFORELP /**< propagation method is called before solving LP */),
#else
				  SCIP_PROPTIMING_BEFORELP, /**< propagation method is called before solving LP */
				  SCIP_PRESOLTIMING_FAST),
#endif
	  model_(NULL),
	  bdsub_(NULL),
	  nvars_(0),
	  vars_(NULL),
	  naux_(0),
	  probability_(NULL)
{
#ifdef BENDERS_PROFILE
	/** initialize statistics */
	names_statistics_.push_back("sepaBenders");
	names_statistics_.push_back("checkBenders");
	for (unsigned i = 0; i < names_statistics_.size(); ++i)
	{
		time_statistics_[names_statistics_[i]] = 0.0;
		count_statistics_[names_statistics_[i]] = 0;
	}
#endif
}

SCIPconshdlrBenders::~SCIPconshdlrBenders()
{
#ifdef BENDERS_PROFILE
	write_statistics();
#endif
}

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
	*result = SCIP_FEASIBLE;
	/**
	 * TODO: DRO needs to pass probability_ to sepaBenders
	 */
	SCIP_CALL(sepaBenders(scip, conshdlr, NULL, result));
	DSPdebugMessage("scip_enfolp: results %d stage %d\n", *result, SCIPgetStage(scip));

	return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for pseudo solutions */
SCIP_DECL_CONSENFOPS(SCIPconshdlrBenders::scip_enfops)
{
	*result = SCIP_FEASIBLE;
	SCIP_CALL(sepaBenders(scip, conshdlr, NULL, result));
	if (*result == SCIP_SEPARATED) *result = SCIP_INFEASIBLE;
	DSPdebugMessage("scip_enfops: results %d stage %d\n", *result, SCIPgetStage(scip));
	return SCIP_OKAY;
}

/** feasibility check method of constraint handler for primal solutions */
SCIP_DECL_CONSCHECK(SCIPconshdlrBenders::scip_check)
{
	*result = SCIP_FEASIBLE;
	SCIP_CALL(checkBenders(scip, conshdlr, sol, result));
	DSPdebugMessage("scip_check: results %d stage(%d)\n", *result, SCIPgetStage(scip));
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

SCIP_DECL_CONSSEPALP(SCIPconshdlrBenders::scip_sepalp)
{
	*result = SCIP_DIDNOTFIND;
	SCIP_CALL(sepaBenders(scip, conshdlr, NULL, result));
	DSPdebugMessage("scip_sepalp: results %d stage %d\n", *result, SCIPgetStage(scip));
	return SCIP_OKAY;
}

SCIP_DECL_CONSSEPASOL(SCIPconshdlrBenders::scip_sepasol)
{
	*result = SCIP_DIDNOTFIND;
	SCIP_CALL(sepaBenders(scip, conshdlr, NULL, result));
	DSPdebugMessage("scip_sepasol: results %d stage %d\n", *result, SCIPgetStage(scip));
	return SCIP_OKAY;
}

SCIP_RETCODE SCIPconshdlrBenders::generate_Benders(
	SCIP *scip,
	SCIP_CONSHDLR *conshdlr,
	SCIP_SOL *sol,
	OsiCuts *cs)
{
	SCIP_Real *vals = NULL; /**< current solution */

	/** allocate memory */
	SCIP_CALL(SCIPallocMemoryArray(scip, &vals, nvars_));

	/** get current solution */
	SCIP_CALL(SCIPgetSolVals(scip, sol, nvars_, vars_, vals));

	/** generate Benders cuts */
	generateCuts(nvars_, vals, cs);

	/** If found Benders cuts */
	for (int i = 0; i < cs->sizeCuts(); ++i)
	{
		/** get cut pointer */
		OsiRowCut *rc = cs->rowCutPtr(i);
		if (!rc)
			continue;

		const CoinPackedVector cutrow = rc->row();
		if (cutrow.getNumElements() == 0)
			rc->setEffectiveness(0.0);
		else
			rc->setEffectiveness(rc->violated(vals) / cutrow.twoNorm());
	}

	/** free memory */
	SCIPfreeMemoryArray(scip, &vals);

	return SCIP_OKAY;
}

SCIP_RETCODE SCIPconshdlrBenders::sepaBenders(
	SCIP *scip,
	SCIP_CONSHDLR *conshdlr,
	SCIP_SOL *sol,
	SCIP_RESULT *result)
{
	double stime = CoinGetTimeOfDay();

	/**< Benders cut placeholder */
	OsiCuts cs;

	/** generate Benders cuts */
	SCIP_CALL(generate_Benders(scip, conshdlr, sol, &cs));

	/** If found Benders cuts */
	for (int i = 0; i < cs.sizeCuts(); ++i)
	{
		/** get cut pointer */
		OsiRowCut * rc = cs.rowCutPtr(i);
		if (!rc) continue;

		const CoinPackedVector cutrow = rc->row();
		if (cutrow.getNumElements() == 0) continue;

#ifdef DSP_DEBUG
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
#endif

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

#ifdef DSP_DEBUG
			DSPdebugMessage("found Benders (%s) cut: act=%f, lhs=%f, norm=%f, eff=%f, min=%f, max=%f (range=%f)\n",
				isOptimalityCut ? "opti" : "feas",
				SCIPgetRowLPActivity(scip, row), SCIProwGetLhs(row), SCIProwGetNorm(row),
				SCIPgetCutEfficacy(scip, sol, row),
				SCIPgetRowMinCoef(scip, row), SCIPgetRowMaxCoef(scip, row),
				SCIPgetRowMaxCoef(scip, row)/SCIPgetRowMinCoef(scip, row));
			DSPdebug(rc->print());
#endif

			/** flush all changes before adding cut */
			SCIP_CALL(SCIPflushRowExtensions(scip, row));

			if (SCIPisGT(scip, rc->effectiveness(), 0.0))
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

			/** add cut to global pool */
			SCIP_CALL(SCIPaddPoolCut(scip, row));
			DSPdebugMessage("number of cuts in global cut pool: %d\n", SCIPgetNPoolCuts(scip));

			/** release the row */
			SCIP_CALL(SCIPreleaseRow(scip, &row));
		}
	}

#ifdef BENDERS_PROFILE
	time_statistics_["sepaBenders"] += CoinGetTimeOfDay() - stime;
	count_statistics_["sepaBenders"]++;
#endif

	return SCIP_OKAY;
}

SCIP_RETCODE SCIPconshdlrBenders::checkBenders(
	SCIP *scip,
	SCIP_CONSHDLR *conshdlr,
	SCIP_SOL *sol,
	SCIP_RESULT *result)
{
	double stime = CoinGetTimeOfDay();

	/**< Benders cut placeholder */
	OsiCuts cs;

	/** generate Benders cuts */
	SCIP_CALL(generate_Benders(scip, conshdlr, sol, &cs));

	/** If found Benders cuts */
	for (int i = 0; i < cs.sizeCuts(); ++i)
	{
		/** get cut pointer */
		OsiRowCut *rc = cs.rowCutPtr(i);
		if (!rc)
			continue;

		const CoinPackedVector cutrow = rc->row();
		if (cutrow.getNumElements() == 0)
			continue;

		if (SCIPisGT(scip, rc->effectiveness(), 0.0))
		{
			*result = SCIP_INFEASIBLE;
			break;
		}
	}

#ifdef BENDERS_PROFILE
	time_statistics_["checkBenders"] += CoinGetTimeOfDay() - stime;
	count_statistics_["checkBenders"]++;
#endif

	return SCIP_OKAY;
}

void SCIPconshdlrBenders::setBdSub(BdSub * bdsub) {
	bdsub_ = bdsub;
	FREE_ARRAY_PTR(probability_);
	probability_ = new double[model_->getNumSubproblems()];

	if (isStochastic()) {
		// extract stochastic model
		TssModel* tss = dynamic_cast<TssModel*>(model_);
		CoinCopyN(tss->getProbability(), tss->getNumScenarios(), probability_);
	} else {
		CoinFillN(probability_, model_->getNumSubproblems(), 1.0);
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

void SCIPconshdlrBenders::generateCuts(
	int size,  /**< [in] size of x */
	double *x, /**< [in] master solution */
	OsiCuts *cuts /**< [out] cuts generated */)
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

void SCIPconshdlrBenders::write_statistics()
{
	fstream fp;
	fp.open("Benders_statistics.txt", ios::out);
	for (unsigned i = 0; i < names_statistics_.size(); ++i)
	{
		fp << names_statistics_[i] << ","
		   << time_statistics_[names_statistics_[i]] << ","
		   << count_statistics_[names_statistics_[i]] << endl;
	}
	fp.close();
}
