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

	/** subproblems */
	TssModel * tssmodel = NULL;

	BGN_TRY_CATCH

	/** allocate memory */
	cutval = new double * [bdsub_->getNumSubprobs()];
	cutrhs = new double [bdsub_->getNumSubprobs()];

	int ncols = model_->getNumCouplingCols();

	/** retrieve tssmodel */
	tssmodel = dynamic_cast<TssModel*>(model_);
	if (!tssmodel) throw "Invalid model type cast";

	/** retrieve probability */
	const double * probability = tssmodel->getProbability();

	/** retrieve dense vector of solution */
	x = solution->denseVector(ncols);

	/** generate cut at the solution */
	bdsub_->generateCuts(ncols, x, cutval, cutrhs);

	for (int j = 0; j < bdsub_->getNumSubprobs(); ++j)
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
		case DSP_STAT_FEASIBLE:
			cuttype.push_back(None);
			break;
		default:
			throw "Unexpected return from Benders cut generation.";
			break;
		}

		DSPdebugMessage("Subproblem index %d probability %e\n", parProcIdx_[j], probability[parProcIdx_[j]]);
		if (bdsub_->getStatus(j) == DSP_STAT_PRIM_INFEASIBLE ||
			bdsub_->getStatus(j) == DSP_STAT_OPTIMAL)
		{
			vec.clear();
			for (int k = 0; k < ncols; ++k)
				if (fabs(cutval[j][k]) > 1.0e-8)
					vec.insert(k, cutval[j][k] * probability[parProcIdx_[j]]);

			OsiRowCut cut;
			cut.setRow(vec);
			cut.setUb(COIN_DBL_MAX);
			cut.setLb(cutrhs[j] * probability[parProcIdx_[j]]);
			DSPdebug(cut.print());
			cuts->insert(cut);
		}
	}

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}
