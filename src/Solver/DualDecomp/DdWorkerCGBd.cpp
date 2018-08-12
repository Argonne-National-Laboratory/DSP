/*
 * DdWorkerCGBd.cpp
 *
 *  Created on: May 23, 2016
 *      Author: kibaekkim
 */

//#define DSP_DEBUG

#include "DdWorkerCGBd.h"

DdWorkerCGBd::DdWorkerCGBd(DspParams * par, DecModel * model, DspMessage * message):
	DdWorkerCG(par, model, message), bdsub_(NULL) {
	// TODO Auto-generated constructor stub
}

DdWorkerCGBd::~DdWorkerCGBd() {
	FREE_PTR(bdsub_);
}

DSP_RTN_CODE DdWorkerCGBd::init() {

	/** subproblems */
	TssModel * tssmodel = NULL;

	BGN_TRY_CATCH

	tssmodel = dynamic_cast<TssModel*>(model_);
	if (!tssmodel) throw "Invalid model type cast";

    DSPdebugMessage("This has %d subproblem indices.\n", parProcIdxSize_);
    DSPdebug(DspMessage::printArray(parProcIdxSize_,parProcIdx_));

	bdsub_ = new BdSub(par_);
	DSP_RTN_CHECK_THROW(bdsub_->setSubIndices(parProcIdxSize_, parProcIdx_));
	DSP_RTN_CHECK_THROW(bdsub_->loadProblem(tssmodel));

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdWorkerCGBd::generateCuts(
	CoinPackedVector * solution, /**< solutions to evaluate */
	OsiCuts *          cuts,     /**< cuts */
	vector<int> &      cuttype   /**< cut type: Opt or Feas */) {
#define FREE_MEMORY                                     \
	FREE_2D_ARRAY_PTR(bdsub_->getNumSubprobs(), cutval) \
	FREE_ARRAY_PTR(cutrhs)                              \
	FREE_ARRAY_PTR(x)

	double ** cutval = NULL;
	double *  cutrhs = NULL;
	double *  x = NULL;

	CoinPackedVector vec;

	BGN_TRY_CATCH

	/** allocate memory */
	cutval = new double * [bdsub_->getNumSubprobs()];
	cutrhs = new double [bdsub_->getNumSubprobs()];

	int ncols = model_->getNumCouplingCols();

	/** retrieve dense vector of solution */
	x = solution->denseVector(ncols);
	DSPdebug(DspMessage::printArray(ncols, x));

	/** generate cut at the solution */
	bdsub_->generateCuts(ncols, x, cutval, cutrhs);

	/***********************************************************
	 * Probability is already applied to the recourse problems.
	 ***********************************************************/

	for (int j = bdsub_->getNumSubprobs() - 1; j >= 0; --j)
	{
		switch (bdsub_->getStatus(j))
		{
		case DSP_STAT_PRIM_INFEASIBLE:
			DSPdebugMessage("Generating feasibility cut.\n");
			cuttype.push_back(Feas);
			break;
		case DSP_STAT_OPTIMAL:
			DSPdebugMessage("Generating optimality cut.\n");
			cuttype.push_back(Opt);
			break;
		default:
			message_->print(0, "Benders subproblem %d returns code %d.\n", j, bdsub_->getStatus(j));
			//throw "Unexpected return from Benders cut generation.";
			break;
		}

		int s = bdsub_->getSubprobIndex(j);
		DSPdebugMessage("Subproblem index %d\n", s);
		if (bdsub_->getStatus(j) == DSP_STAT_PRIM_INFEASIBLE ||
			bdsub_->getStatus(j) == DSP_STAT_OPTIMAL)
		{
			vec.clear();
			for (int k = 0; k < ncols; ++k)
				if (fabs(cutval[j][k]) > 1.0e-8)
					vec.insert(k, cutval[j][k]);

			OsiRowCut cut;
			cut.setRow(vec);
			cut.setUb(COIN_DBL_MAX);
			cut.setLb(cutrhs[j]);
			DSPdebug(cut.print());
			cuts->insert(cut);
		}

		if (bdsub_->getStatus(j) == DSP_STAT_PRIM_INFEASIBLE)
			break;
	}

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}
