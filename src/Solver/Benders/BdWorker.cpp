/*
 * BdWorker.cpp
 *
 *  Created on: Feb 18, 2016
 *      Author: kibaekkim
 */

//#define DSP_DEBUG
#include "Solver/Benders/BdWorker.h"
#include "CoinHelperFunctions.hpp"

BdWorker::BdWorker(
		DecModel *   model,  /**< model pointer */
		DspParams *  par,    /**< parameters */
		DspMessage * message /**< message pointer */) :
model_(model),
par_(par),
message_(message) {

	/** parameters */
	parProcIdxSize_ = par_->getIntPtrParamSize("ARR_PROC_IDX");
	parProcIdx_     = par_->getIntPtrParam("ARR_PROC_IDX");

	/** create Benders subproblem solver */
	bdsub_ = new BdSub(par_);
	bdsub_->setSubIndices(parProcIdxSize_, parProcIdx_);
	bdsub_->loadProblem(model_);
}

BdWorker::BdWorker(const BdWorker& rhs) :
model_(rhs.model_),
par_(rhs.par_),
message_(rhs.message_),
parProcIdxSize_(rhs.parProcIdxSize_),
parProcIdx_(rhs.parProcIdx_) {
	bdsub_ = rhs.bdsub_->clone();
}

BdWorker::~BdWorker() {
	FREE_PTR(bdsub_);
}

DSP_RTN_CODE BdWorker::generateCuts(int nx, int naux, const double* x, OsiCuts& cs) {

#define FREE_MEMORY \
	FREE_2D_ARRAY_PTR(parProcIdxSize_,cutval) \
	FREE_ARRAY_PTR(cutrhs)

	double** cutval = NULL;
	double*  cutrhs = NULL;

	BGN_TRY_CATCH

	/** allocate memory for cuts */
	cutval = new double * [parProcIdxSize_];
	cutrhs = new double [parProcIdxSize_];
	for (int s = 0; s < parProcIdxSize_; ++s)
		cutval[s] = NULL;
	CoinZeroN(cutrhs, parProcIdxSize_);

	/** generate cuts */
	bdsub_->generateCuts(nx, x, cutval, cutrhs);

	/** collect cuts */
	collectCuts(nx, naux, cutval, cutrhs, cs);

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

DSP_RTN_CODE BdWorker::collectCuts(
		int nx,
		int naux,
		double** cut,
		double* rhs,
		OsiCuts& cs) {

#define FREE_MEMORY \
	FREE_2D_ARRAY_PTR(naux,aggcut) \
	FREE_ARRAY_PTR(aggrhs)

	double** aggcut = NULL;
	double* aggrhs = NULL;
	CoinPackedVector vec;

	BGN_TRY_CATCH

	/** is there a feasibility cut? */
	int fcut = -1;
	for (int i = 0; i < bdsub_->getNumSubprobs(); ++i)
		if (bdsub_->getStatus(i) == DSP_STAT_PRIM_INFEASIBLE) {
			fcut = i;
			break;
		}

	if (fcut > -1) {
		/** initialize vector */
		vec.clear();
		/** set it as sparse */
		for (int j = 0; j < nx; ++j) {
			if (fabs(cut[fcut][j]) > 1e-10)
				vec.insert(j, cut[fcut][j]);
		}

		/** create row cut */
		OsiRowCut rc;
		rc.setRow(vec);
		rc.setUb(COIN_DBL_MAX);
		rc.setLb(rhs[fcut]);

		/** store cut */
		cs.insert(rc);
	} else {
		/** total number of subproblems */
		int nsubprobs = parProcIdxSize_;

		/** allocate memory for cuts */
		aggcut = new double * [naux];
		aggrhs = new double [naux];
		for (int s = 0; s < naux; ++s) {
			aggcut[s] = new double [nx];
			CoinZeroN(aggcut, nx);
		}
		CoinZeroN(aggrhs, naux);

		/** aggregate cuts */
		for (int s = 0; s < nsubprobs; ++s) {
			int i = s % naux;
			for (int j = 0; j < nx; ++j)
				aggcut[i][j] += cut[s][j];
			aggrhs[i] += rhs[s];
		}

		/** construct cuts */
		for (int i = 0; i < naux; ++i) {

			/** initialize vector */
			vec.clear();

			/** set it as sparse */
			for (int j = 0; j < nx; ++j) {
				if (fabs(aggcut[i][j]) > 1e-10)
					vec.insert(j, aggcut[i][j]);
			}
			vec.insert(nx + i, 1.0);

			/** create row cut */
			OsiRowCut rc;
			rc.setRow(vec);
			rc.setUb(COIN_DBL_MAX);
			rc.setLb(aggrhs[i]);

			/** store cut */
			cs.insert(rc);
		}
	}

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}
