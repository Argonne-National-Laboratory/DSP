/*
 * TssEval.cpp
 *
 *  Created on: Mar 16, 2015
 *      Author: kibaekkim
 */

#define DSP_DEBUG

#include "Solver/TssEval.h"
#include "Solver/TssBdSub.h"
#include <Utility/DspMessage.h>

TssEval::TssEval() :
hasSolution_(false),
parNumCores_(1)
{
	/** nothing to do */
}

TssEval::~TssEval()
{
	/** nothing to do */
}

/** solve */
DSP_RTN_CODE TssEval::solve()
{
#define FREE_MEMORY                                              \
	FREE_ARRAY_PTR(objval_reco);                                 \
	FREE_2D_ARRAY_PTR(model_->getNumScenarios(), solution_reco); \
	FREE_PTR(sub);

	double * objval_reco = NULL;
	double ** solution_reco = NULL;
	TssBdSub * sub = NULL;

	BGN_TRY_CATCH

	/** parameters */
	parNumCores_ = par_->getIntParam("BD/NUM_CORES");

	/** allocate memory */
	objval_reco = new double [model_->getNumScenarios()];
	solution_reco = new double * [model_->getNumScenarios()];
	for (int s = 0; s < model_->getNumScenarios(); ++s)
		solution_reco[s] = new double [model_->getNumCols(1)];
	sub = new TssBdSub(par_);
	for (int s = 0; s < model_->getNumScenarios(); ++s)
		sub->scenarios_.push_back(s);

	/** load subproblems */
	sub->loadProblem(model_, 0, NULL, 0);

	/** evaluate solution */
#ifdef DSP_DEBUG2
	for (int j = 0; j < model_->getNumCols(0); ++j)
		printf("x%d %e\n", j, solution_[j]);
#endif
	sub->solveRecourse(solution_, objval_reco, solution_reco, parNumCores_);

	/** get objective value */
	primalBound_ = 0.0;
	for (int j = 0; j < model_->getNumCols(0); ++j)
		primalBound_ += model_->getObjCore(0)[j] * solution_[j];
	for (int s = 0; s < model_->getNumScenarios(); ++s)
	{
		DSPdebugMessage("scenario %d objective %e\n", s, objval_reco[s]);
		primalBound_ += objval_reco[s] * model_->getProbability()[s];
	}
	DSPdebugMessage("primalBound %e\n", primalBound_);

	/** collect second-stage solutions */
	for (int s = 0; s < model_->getNumScenarios(); ++s)
	{
		CoinCopyN(solution_reco[s], model_->getNumCols(1),
				solution_ + model_->getNumCols(0) + model_->getNumCols(1) * s);
	}

	END_TRY_CATCH_RTN(FREE_MEMORY,STO_RTN_ERR);

	/** free memory */
	FREE_MEMORY;

	return STO_RTN_OK;
#undef FREE_MEMORY
}

/** set solution to evaluate */
void TssEval::setSolution(double * solution)
{
	if (solution == NULL) return;

	hasSolution_ = true;
	CoinCopyN(solution, model_->getNumCols(0), solution_);
}

