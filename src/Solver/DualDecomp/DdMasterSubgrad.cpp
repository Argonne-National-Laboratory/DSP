/*
 * DdMasterSubgrad.cpp
 *
 *  Created on: Feb 16, 2016
 *      Author: kibaekkim
 */

#include "Model/DecTssModel.h"
#include "Solver/DualDecomp/DdMasterSubgrad.h"

DdMasterSubgrad::DdMasterSubgrad(
		DspParams *  par,     /**< parameter pointer */
		DecModel *   model,   /**< model pointer */
		DspMessage * message, /**< message pointer */
		int nworkers          /**< number of workers */):
DdMasterSync(par, model, message, nworkers),
nstalls_(0),
stepscal_(2.0),
stepsize_(0.0),
gradient_(NULL),
multipliers_(NULL) {}

DdMasterSubgrad::~DdMasterSubgrad()
{
	FREE_ARRAY_PTR(gradient_);
	FREE_ARRAY_PTR(multipliers_);
}

DSP_RTN_CODE DdMasterSubgrad::init()
{
	BGN_TRY_CATCH

	DdMasterSync::init();

	status_ = DSP_STAT_OPTIMAL;

	/** create problem */
	createProblem();

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMasterSubgrad::solve()
{
	BGN_TRY_CATCH

	double cputime  = CoinCpuTime();
	double walltime = CoinGetTimeOfDay();

	/** get new Lagrangian multipliers */
	for (int j = 0; j < model_->getNumCouplingRows(); ++j)
	{
		multipliers_[j] += stepsize_ * gradient_[j];
		if (model_->getSenseCouplingRow(j) == 'L') /** <= : nonnegative multipliers */
			multipliers_[j] = CoinMax((double) 0.0, multipliers_[j]);
		else if (model_->getSenseCouplingRow(j) == 'G') /** >= : nonpositive multipliers */
			multipliers_[j] = CoinMin((double) 0.0, multipliers_[j]);
	}

	/** store final multipliers */
	DecTssModel * decTssModel = dynamic_cast<DecTssModel*>(model_);
	if (model_->nonanticipativity() && decTssModel != NULL)
	{
		/** nonanticipativity constraints must be converted back from a different representation */
		decTssModel->convertLagrangianFromAlternative(multipliers_, primsol_);
	}
	else
	{
		/** copy Lagrangian multipliers */
		for (int j = 0; j < model_->getNumCouplingRows(); ++j)
			primsol_[j] = multipliers_[j];
	}

	/** update statistics */
	double * s_primsol = new double [model_->getNumCouplingRows()];
	CoinCopyN(primsol_, model_->getNumCouplingRows(), s_primsol);
	s_primsols_.push_back(s_primsol);
	s_primsol = NULL;
	s_cputimes_.push_back(CoinCpuTime() - cputime);
	s_walltimes_.push_back(CoinGetTimeOfDay() - walltime);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMasterSubgrad::createProblem()
{
	BGN_TRY_CATCH

	/** allocate memory */
	primsol_     = new double [model_->getNumCouplingRows()];
	gradient_    = new double [model_->getNumCouplingRows()];
	multipliers_ = new double [model_->getNumCouplingRows()];

	/** initialize values */
	CoinZeroN(primsol_, model_->getNumCouplingRows());
	CoinZeroN(gradient_, model_->getNumCouplingRows());
	CoinZeroN(multipliers_, model_->getNumCouplingRows());

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMasterSubgrad::updateProblem()
{
	BGN_TRY_CATCH

	/** compute gradient */
	DecTssModel * decTssModel = dynamic_cast<DecTssModel*>(model_);
	if (model_->nonanticipativity() && decTssModel != NULL)
	{
		/** nonanticipativity constraints must be converted to a different representation */
		for (int i = 0; i < model_->getNumCouplingRows(); i++)
			gradient_[i] = decTssModel->evalLhsCouplingRowAlternative(i, subsolution_);
	}
	else
	{
		for (int i = 0; i < model_->getNumCouplingRows(); i++)
			gradient_[i] = model_->evalLhsCouplingRow(i, subsolution_) - model_->getRhsCouplingRow(i);
	}

	/** Lagrangian value */
	double newobj = 0.0;
	for (int s = 0; s < nsubprobs_; ++s)
		newobj += subdualobj_[s];

	/** update dual bound */
	DSPdebugMessage("-> Lagrangian value %e\n", newobj);
	if (newobj > bestdualobj_)
	{
		bestdualobj_ = newobj;
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
	for (int j = 0; j < model_->getNumCouplingRows(); ++j)
		denom += gradient_[j] * gradient_[j];
	if (bestprimobj_ < 1.0e+20)
		stepsize_ = stepscal_ * (bestprimobj_ - newobj) / denom;
	else
		stepsize_ = stepscal_;
	DSPdebugMessage("-> step size %e\n", stepsize_);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}
