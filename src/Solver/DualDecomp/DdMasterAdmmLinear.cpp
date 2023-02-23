/*
 * DdMasterAdmmLinear.cpp
 *
 *  Created on: Feb 10, 2023
 *      Author: hideakiv
 */

#define DSP_DEBUG

#include "Model/DecTssModel.h"
#include "Solver/DualDecomp/DdMasterAdmmLinear.h"

DdMasterAdmmLinear::DdMasterAdmmLinear(
		DecModel *   model,  /**< model pointer */
		DspParams *  par,    /**< parameter pointer */
		DspMessage * message /**< message pointer */):
DdMaster(model, par, message),
itercount_(0),
nstalls_(0),
stepscal_(2.0),
stepsize_(0.0),
gradient_(NULL),
multipliers_(NULL),
isPrimFeas(true),
prev_gradient_(NULL),
mean_multipliers_(NULL) {}

DdMasterAdmmLinear::DdMasterAdmmLinear(const DdMasterAdmmLinear& rhs) :
DdMaster(rhs),
itercount_(rhs.itercount_),
nstalls_(rhs.nstalls_),
stepscal_(rhs.stepscal_),
stepsize_(rhs.stepsize_),
isPrimFeas(rhs.isPrimFeas) {
	gradient_         = new double [model_->getNumCouplingRows()];
	multipliers_      = new double [model_->getNumCouplingRows()];
    DecTssModel * decTssModel = dynamic_cast<DecTssModel*>(model_);
    prev_gradient_    = new double [decTssModel->getNumCols(0)];
    mean_multipliers_ = new double [decTssModel->getNumCols(0)]; 
	CoinCopyN(rhs.gradient_, model_->getNumCouplingRows(), gradient_);
	CoinCopyN(rhs.multipliers_, model_->getNumCouplingRows(), multipliers_);
	CoinCopyN(rhs.prev_gradient_, decTssModel->getNumCols(0), prev_gradient_);
    CoinCopyN(rhs.mean_multipliers_, decTssModel->getNumCols(0), mean_multipliers_);
}

DdMasterAdmmLinear::~DdMasterAdmmLinear()
{
	FREE_ARRAY_PTR(gradient_);
	FREE_ARRAY_PTR(multipliers_);
	FREE_ARRAY_PTR(prev_gradient_);
    FREE_ARRAY_PTR(mean_multipliers_);
}

DSP_RTN_CODE DdMasterAdmmLinear::init()
{
	BGN_TRY_CATCH

	DdMaster::init();

	status_ = DSP_STAT_OPTIMAL;

	/** create problem */
	createProblem();

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMasterAdmmLinear::solve()
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
    //
    assert(model_->isStochastic());
    //
	if (model_->isStochastic() && decTssModel != NULL)
	{
		/** copy Lagrangian multipliers */
		CoinCopyN(multipliers_, model_->getNumCouplingRows(), &primsol_[0]);
	}
	else
	{
        /* TODO: handle non-stochastic cases */
        assert(false);
		/** copy Lagrangian multipliers */
		CoinCopyN(multipliers_, model_->getNumCouplingRows(), &primsol_[0]);
	}

    assert(decTssModel != NULL);
    /** take average of lagrange multipliers */
    for (int j = 0; j < decTssModel->getNumCols(0); j++) {
        mean_multipliers_[j] = 0.0;
        for (int s = 0; s < model_->getNumSubproblems(); s++) {
            int i = s * decTssModel->getNumCols(0) + j;
            mean_multipliers_[j] += primsol_[i];
        }
        mean_multipliers_[j] /= model_->getNumSubproblems();
	    DSPdebugMessage("-> Mean multipliers %e\n", mean_multipliers_[j]);
    }
    for (int j = 0; j < decTssModel->getNumCols(0); j++) {
        if (mean_multipliers_[j] > 1e-10) {
            isPrimFeas = false;
            break;
        }
    }

    /** update lagrange multipliers */
    for (int i = 0; i < model_->getNumCouplingRows(); ++i)
	{
        int j = i % decTssModel->getNumCols(0);
		primsol_[i] += mean_multipliers_[j] - stepsize_ * prev_gradient_[j];
	}

    /** store average gradients */
    for (int j = 0; j < decTssModel->getNumCols(0); j++) {
        prev_gradient_[j] = 0.0;
        for (int s = 0; s < model_->getNumSubproblems(); s++) {
            int i = s * decTssModel->getNumCols(0) + j;
            prev_gradient_[j] += gradient_[i];
        }
        prev_gradient_[j] /= model_->getNumSubproblems();
    }

	/** retrieve lambda */
	lambda_ = primsol_;

	/** update statistics */
	double * s_primsol = new double [model_->getNumCouplingRows()];
	CoinCopyN(&primsol_[0], model_->getNumCouplingRows(), s_primsol);
	s_primsols_.push_back(s_primsol);
	s_primsol = NULL;
	s_cputimes_.push_back(CoinCpuTime() - cputime);
	s_walltimes_.push_back(CoinGetTimeOfDay() - walltime);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMasterAdmmLinear::createProblem()
{
	BGN_TRY_CATCH
	DecTssModel * decTssModel = dynamic_cast<DecTssModel*>(model_);

	/** allocate memory */
	gradient_         = new double [model_->getNumCouplingRows()];
	multipliers_      = new double [model_->getNumCouplingRows()];
	prev_gradient_    = new double [decTssModel->getNumCols(0)];
    mean_multipliers_ = new double [decTssModel->getNumCols(0)];

	/** initialize values */
	primsol_.resize(model_->getNumCouplingRows());
	CoinZeroN(gradient_, model_->getNumCouplingRows());
	CoinZeroN(multipliers_, model_->getNumCouplingRows());
	CoinZeroN(prev_gradient_, decTssModel->getNumCols(0));
	CoinZeroN(mean_multipliers_, decTssModel->getNumCols(0));

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMasterAdmmLinear::updateProblem()
{
	BGN_TRY_CATCH

	/** compute gradient */
	DecTssModel * decTssModel = dynamic_cast<DecTssModel*>(model_);
    //
    assert(model_->isStochastic());
    //
	if (model_->isStochastic() && decTssModel != NULL)
	{
		for (int i = 0; i < model_->getNumCouplingRows(); i++) {
            int s = i / decTssModel->getNumCols(0);
	        int j = i % decTssModel->getNumCols(0);
	        gradient_[i] = subsolution_[s][j];
        }
	}
	else
	{
        /* TODO: handle non-stochastic cases */
        assert(false);
		for (int i = 0; i < model_->getNumCouplingRows(); i++)
			gradient_[i] = model_->evalLhsCouplingRow(i, subsolution_) - model_->getRhsCouplingRow(i);
	}

	/** Augmented Lagrangian value */

	double newobj = 0.0;
    if (isPrimFeas) {
        for (int s = 0; s < model_->getNumSubproblems(); ++s)
		    newobj += subdualobj_[s];
    } else {
        newobj = -COIN_DBL_MAX;
    }

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
    itercount_++;

	/** calculate step size */
	switch (par_->getIntParam("DD/MASTER_STEP_RULE"))
    {
	case Polyak:
		double denom;
        denom = 0.0;
		for (int j = 0; j < model_->getNumCouplingRows(); ++j)
			denom += gradient_[j] * gradient_[j];
		if (bestprimobj_ < 1.0e+20 && newobj > -1.0e+20)
			stepsize_ = stepscal_ * (bestprimobj_ - newobj) / denom;
		else
			stepsize_ = stepscal_;
		break;
	case SSNS:
		stepsize_ = 1.0 / ( 0.0 + itercount_);
        break;
	}
	
	DSPdebugMessage("-> step size %e\n", stepsize_);

    /** update statistics */
    s_statuses_.push_back(DSP_STAT_OPTIMAL);
	s_primobjs_.push_back(bestprimobj_);
	s_dualobjs_.push_back(newobj);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}
