/*
 * TssDdMasterSubgrad.cpp
 *
 *  Created on: Mar 3, 2015
 *      Author: kibaekkim
 */

//#define DSP_DEBUG

#include "Solver/TssDdMasterSubgrad.h"
#include "Utility/StoMessage.h"

TssDdMasterSubgrad::~TssDdMasterSubgrad()
{
	FREE_ARRAY_PTR(solution_);
	FREE_ARRAY_PTR(gradient_);
	FREE_ARRAY_PTR(multipliers_);
	ncols_ = 0;
}

/** create problem */
STO_RTN_CODE TssDdMasterSubgrad::createProblem(const TssModel * model)
{
	nscenarios_ = model->getNumScenarios();
	ncols_first_ = model->getNumCols(0);
	ncols_ = nscenarios_ * ncols_first_;

	/** allocate memory */
	solution_ = new double [ncols_];
	gradient_ = new double [ncols_];
	multipliers_ = new double [ncols_];

	/** initialize values */
	CoinZeroN(solution_, ncols_);
	CoinZeroN(gradient_, ncols_);
	CoinZeroN(multipliers_, ncols_);

	return STO_RTN_OK;
}

/** update problem: may update dual bound */
STO_RTN_CODE TssDdMasterSubgrad::updateProblem(
		double primal_bound, /**< primal bound of the original problem */
		double & dual_bound, /**< dual bound of the original problem */
		double * objvals,    /**< objective values of subproblems */
		double ** solution   /**< subproblem solutions */)
{
	/** copy gradient */
	for (int s = 0; s < nscenarios_; ++s)
		for (int j = 0; j < ncols_first_; ++j)
		{
			if (j + 1 < ncols_first_)
				gradient_[s * ncols_first_ + j] = solution[s][j] - solution[s][j+1];
			else
				gradient_[s * ncols_first_ + j] = solution[s][j] - solution[s][0];
		}

	/** Lagrangian value */
	double newobj = 0.0;
	for (int s = 0; s < nscenarios_; ++s)
		newobj += objvals[s];

	/** update dual bound */
	DSPdebugMessage("-> Lagrangian value %e\n", newobj);
	if (newobj > dual_bound)
	{
		dual_bound = newobj;
		DSPdebugMessage("-> Updated dual bound %e\n", newobj);
		nstalls_ = 0;
	}
	else
	{
		if (nstalls_ < 5)
			nstalls_++;
		else
		{
			stepscal_ *= 0.5;
			DSPdebugMessage("Updated constant parameter %e\n", stepscal_);
			nstalls_ = 0;
		}
	}

	/** calculate step size */
	double denom = 0.0;
	for (int j = 0; j < ncols_; ++j)
		denom += gradient_[j] * gradient_[j];
	if (primal_bound < 1.0e+20)
		stepsize_ = stepscal_ * (primal_bound - newobj) / denom;
	else
		stepsize_ = 1.0;
	DSPdebugMessage("-> step size %e\n", stepsize_);

	return STO_RTN_OK;
}

/** solve problem */
STO_RTN_CODE TssDdMasterSubgrad::solve()
{
	/** get new Lagrangian multipliers */
	for (int j = 0; j < ncols_; ++j)
		multipliers_[j] += stepsize_ * gradient_[j];

	for (int s = 0; s < nscenarios_; ++s)
		for (int j = 0; j < ncols_first_; ++j)
		{
			if (s == 0)
				solution_[s * ncols_first_ + j] = multipliers_[s * ncols_first_ + j] - multipliers_[(nscenarios_ - 1) * ncols_first_ + j];
			else
				solution_[s * ncols_first_ + j] = multipliers_[s * ncols_first_ + j] - multipliers_[(s - 1) * ncols_first_ + j];
		}

	return STO_RTN_OK;
}

