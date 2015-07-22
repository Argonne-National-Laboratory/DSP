/*
 * DecDdMasterSubgrad.cpp
 *
 *  Created on: Mar 3, 2015
 *      Author: kibaekkim, ctjandra
 */

//#define DSP_DEBUG

#include "Solver/DecDdMasterSubgrad.h"
#include "Utility/StoMessage.h"
#include "Model/DecTssModel.h"

DecDdMasterSubgrad::~DecDdMasterSubgrad()
{
	FREE_ARRAY_PTR(solution_);
	FREE_ARRAY_PTR(gradient_);
	FREE_ARRAY_PTR(multipliers_);
	ncoupling_ = 0;
}

/** create problem */
STO_RTN_CODE DecDdMasterSubgrad::createProblem(DecModel * model)
{
	model_ = model;
	nsubprobs_ = model->getNumSubproblems();
	ncoupling_ = model->getNumCouplingRows();
	nonanticipativity_ = model->nonanticipativity();

	/** allocate memory */
	solution_ = new double [ncoupling_];
	gradient_ = new double [ncoupling_];
	multipliers_ = new double [ncoupling_];

	/** initialize values */
	CoinZeroN(solution_, ncoupling_);
	CoinZeroN(gradient_, ncoupling_);
	CoinZeroN(multipliers_, ncoupling_);

	return STO_RTN_OK;
}

/** update problem: may update dual bound */
STO_RTN_CODE DecDdMasterSubgrad::updateProblem(
		double primal_bound, /**< primal bound of the original problem */
		double & dual_bound, /**< dual bound of the original problem */
		double * objvals,    /**< objective values of subproblems */
		double ** solution   /**< subproblem solutions */)
{
	/** compute gradient */
	DecTssModel * decTssModel = dynamic_cast<DecTssModel*>(model_);
	if (nonanticipativity_ && decTssModel != NULL)
	{
		/** nonanticipativity constraints must be converted to a different representation */
		for (int i = 0; i < ncoupling_; i++)
			gradient_[i] = decTssModel->evalLhsRowAlternative(i, solution);
	}
	else
	{
		for (int i = 0; i < ncoupling_; i++)
			gradient_[i] = model_->evalLhsRow(i, solution);
	}

	/** Lagrangian value */
	double newobj = 0.0;
	for (int s = 0; s < nsubprobs_; ++s)
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
	for (int j = 0; j < ncoupling_; ++j)
		denom += gradient_[j] * gradient_[j];
	if (primal_bound < 1.0e+20)
		stepsize_ = stepscal_ * (primal_bound - newobj) / denom;
	else
		stepsize_ = 1.0;
	DSPdebugMessage("-> step size %e\n", stepsize_);

	return STO_RTN_OK;
}

/** solve problem */
STO_RTN_CODE DecDdMasterSubgrad::solve()
{
	/** get new Lagrangian multipliers */
	for (int j = 0; j < ncoupling_; ++j)
		multipliers_[j] += stepsize_ * gradient_[j];

	/** store final multipliers */
	DecTssModel * decTssModel = dynamic_cast<DecTssModel*>(model_);
	if (nonanticipativity_ && decTssModel != NULL)
	{
		/** nonanticipativity constraints must be converted back from a different representation */
		decTssModel->convertLagrangianFromAlternative(multipliers_, solution_);
	}
	else
	{
		/** copy Lagrangian multipliers */
		for (int j = 0; j < ncoupling_; ++j)
			solution_[j] = multipliers_[j];
	}

	return STO_RTN_OK;
}

