#include "Solver/Benders/BendersCallback.h"
#include "Utility/DspMessage.h"
#include "Model/TssModel.h"
#include "OsiCuts.hpp"


/** default constructor */
BendersCallback():
osi_(NULL,)
model_(NULL),
bdsub_(NULL),
nvars_(0),
naux_(0),
probability_(NULL){}

~BendersCallback::BendersCallback(){
	/** empty*/
}

int static BendersCut(void *cbdata, int cbwhere){
	//BGN_TRY_CATCH
	double stime = CoinGetTimeOfDay();

	/**< Benders cut placeholder */
	OsiCuts cs;

	/** generate Benders cuts */
	generate_Benders(osi, &cs);

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

	
#ifdef DSP_DEBUG
			/*DSPdebugMessage("found Benders (%s) cut: act=%f, lhs=%f, norm=%f, eff=%f, min=%f, max=%f (range=%f)\n",
				isOptimalityCut ? "opti" : "feas",
				SCIPgetRowLPActivity(scip, row), SCIProwGetLhs(row), SCIProwGetNorm(row),
				SCIPgetCutEfficacy(scip, sol, row),
				SCIPgetRowMinCoef(scip, row), SCIPgetRowMaxCoef(scip, row),
				SCIPgetRowMaxCoef(scip, row)/SCIPgetRowMinCoef(scip, row));*/
			DSPdebug(rc->print());
#endif

			if (rc->effectiveness()>0.0)
			{
				/** add cut */
				//if (cbwhere==CB_MIPSOL):
					osi->CallbackLazycut(cbdata, rc);
			}
	}

	return DSP_RTN_OK;
}

DSP_RTN_CODE BendersCallback::addBenderscut(DspOsi *osi){

	osi->setCallbackFunc(Benderscut);

}

void BendersCallback::setBdSub(BdSub * bdsub) {
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

DSP_RTN_CODE BendersCallback::generate_Benders(DspOsi * osi, OsiCuts *cs){

	//BGN_TRY_CATCH

	double *vals;
	vals=new double[nvars_];

	/** get current solution */
	osi->cbget(CB_MIPSOL, CB_MIPSOL_SOL, &vals);

	/** generate Benders cuts */
	generateCuts(nvars_, val, cs);

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
	delete [] val;

	return DSP_RTN_OK;
}

DSP_RTN_CODE BendersCallback::sepaBenders(
	DspOsi *osi)
{
	//BGN_TRY_CATCH
	double stime = CoinGetTimeOfDay();

	/**< Benders cut placeholder */
	OsiCuts cs;

	/** generate Benders cuts */
	generate_Benders(osi, &cs));

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

	
#ifdef DSP_DEBUG
			/*DSPdebugMessage("found Benders (%s) cut: act=%f, lhs=%f, norm=%f, eff=%f, min=%f, max=%f (range=%f)\n",
				isOptimalityCut ? "opti" : "feas",
				SCIPgetRowLPActivity(scip, row), SCIProwGetLhs(row), SCIProwGetNorm(row),
				SCIPgetCutEfficacy(scip, sol, row),
				SCIPgetRowMinCoef(scip, row), SCIPgetRowMaxCoef(scip, row),
				SCIPgetRowMaxCoef(scip, row)/SCIPgetRowMinCoef(scip, row));*/
			DSPdebug(rc->print());
#endif

			if (rc->effectiveness()>0.0))
			{
				/** add cut */
				osi->CallbackLazyCut(cbdata, rc);
			}
	}

#ifdef BENDERS_PROFILE
	time_statistics_["sepaBenders"] += CoinGetTimeOfDay() - stime;
	count_statistics_["sepaBenders"]++;
#endif

	return DSP_RTN_OK;
}

 DSP_RTN_CODE BendersCallback::setOriginalVariables(
			int nvars,        /**< number of original variables, including auxiliary variables */
			int         naux  /**< number of auxiliary variables */)
{
	//BGN_TRY_CATCH
	nvars_ = nvars;
	naux_ = naux;

	return DSP_RTN_OK;
}

void BendersCallback::generateCuts(
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

void BendersCallback::aggregateCuts(
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

	//END_TRY_CATCH(FREE_MEMORY)

	FREE_MEMORY

#undef FREE_MEMORY
}
