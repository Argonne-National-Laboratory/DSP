/*
 * DdMasterAdmmLinear.cpp
 *
 *  Created on: Feb 10, 2023
 *      Author: hideakiv
 */

// #define DSP_DEBUG

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
mean_multipliers_(NULL),
auglagrange(0.0),
primres(0.0) {}

DdMasterAdmmLinear::DdMasterAdmmLinear(const DdMasterAdmmLinear& rhs) :
DdMaster(rhs),
itercount_(rhs.itercount_),
nstalls_(rhs.nstalls_),
stepscal_(rhs.stepscal_),
stepsize_(rhs.stepsize_),
isPrimFeas(rhs.isPrimFeas),
auglagrange(rhs.auglagrange),
primres(rhs.primres) {
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
	DecTssModel * decTssModel = dynamic_cast<DecTssModel*>(model_);
    assert(decTssModel != NULL);
    /**     take average of current lagrange multipliers */
    for (int j = 0; j < decTssModel->getNumCols(0); j++) {
        double mean_multiplier = 0.0;
        for (int s = 0; s < model_->getNumSubproblems(); s++) {
            int i = s * decTssModel->getNumCols(0) + j;
            mean_multiplier += multipliers_[i];
        }
        mean_multipliers_[j] = mean_multiplier / model_->getNumSubproblems();
    }
    /** update augmented lagrangian */
    auglagrange = 0.0;
    for (int i = 0; i < model_->getNumCouplingRows(); ++i)
	{
        int j = i % decTssModel->getNumCols(0);
		auglagrange += 0.5 * stepsize_ * (gradient_[i] * gradient_[i] - prev_gradient_[j] * prev_gradient_[j]);
	} // continue to updateProblem()

    /**     update lagrange multipliers */
	for (int i = 0; i < model_->getNumCouplingRows(); ++i)
	{
        int j = i % decTssModel->getNumCols(0);
		multipliers_[i] += stepsize_ * (gradient_[i] - prev_gradient_[j]) - mean_multipliers_[j];
		if (model_->getSenseCouplingRow(i) == 'L') /** <= : nonnegative multipliers */
			multipliers_[i] = CoinMax((double) 0.0, multipliers_[i]);
		else if (model_->getSenseCouplingRow(i) == 'G') /** >= : nonpositive multipliers */
			multipliers_[i] = CoinMin((double) 0.0, multipliers_[i]);
	}

	/** store final multipliers */
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

    /** store average gradients */
    primres = 0.0;
    for (int j = 0; j < decTssModel->getNumCols(0); j++) {
        double mean_grad = 0.0;
        for (int s = 0; s < model_->getNumSubproblems(); s++) {
            int i = s * decTssModel->getNumCols(0) + j;
            mean_grad += gradient_[i];
        }
        mean_grad /= model_->getNumSubproblems();
        
        double prim_res_val = stepsize_ * (prev_gradient_[j] - mean_grad);
	    DSPdebugMessage("-> Updated mean multipliers %e\n", prim_res_val);

        primres += prim_res_val * prim_res_val;
        prev_gradient_[j] = mean_grad;
    }
    primres /= decTssModel->getNumCols(0);
	DSPdebugMessage("-> Primal residual %e\n", primres);

    isPrimFeas = primres < 1e-10;

	/** retrieve lambda */
	lambda_.resize(model_->getNumCouplingRows());
	
	if (model_->isStochastic() && decTssModel != NULL)
	{
		for (int i = 0; i < model_->getNumCouplingCols(); ++i) {
			for (int j = 0; j < model_->getNumSubproblems(); ++j) {
				int k = model_->getNumCouplingCols()*j+i;
				lambda_[k] = primsol_[k] + model_->getCouplingColsObjs()[i] * decTssModel->getProbability()[j];
			}
		}
	} else {
        /* TODO: handle non-stochastic cases */
        assert(false);
		CoinCopyN(&primsol_[0], model_->getNumCouplingRows(), &lambda_[0]);
	}

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
    for (int s = 0; s < model_->getNumSubproblems(); ++s)
		newobj += subdualobj_[s];

    /** update augmented lagrangian */
    auglagrange += newobj;


    if (!isPrimFeas) {
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
            switch (par_->getIntParam("DD/MASTER_STEP_RULE"))
			{
			case Polyak:
				stepscal_ *= 0.5;
			case SSNS:
				stepscal_ *= 1.0;
			}
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
    if (isPrimFeas){
        s_statuses_.push_back(DSP_STAT_OPTIMAL);
    } else {
        s_statuses_.push_back(DSP_STAT_PRIM_INFEASIBLE);
    }
	s_primobjs_.push_back(bestprimobj_);
	s_dualobjs_.push_back(newobj);
    s_auglagrange_.push_back(auglagrange);
    s_primres_.push_back(primres);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

/** write output to a file */
void DdMasterAdmmLinear::write(const char * filename)
{
	BGN_TRY_CATCH

	ofstream myfile;
	myfile.open(filename);
	myfile << "Iter";
	myfile << ",Status";
	myfile << ",Prim";
	myfile << ",Dual";
	myfile << ",Cpu";
	myfile << ",Wall";
    myfile << ",AugLag";
    myfile << ",PrimRes";
	myfile << "\n";
	for (unsigned i = 0; i < s_statuses_.size(); ++i)
	{
		myfile << i;
		myfile << "," << s_statuses_[i];
		myfile << "," << scientific << setprecision(5) << s_primobjs_[i];
		myfile << "," << scientific << setprecision(5) << s_dualobjs_[i];
		myfile << "," << fixed << setprecision(2) << s_cputimes_[i];
		myfile << "," << fixed << setprecision(2) << s_walltimes_[i];
        myfile << "," << fixed << setprecision(2) << s_auglagrange_[i];
        myfile << "," << fixed << setprecision(2) << s_primres_[i];
		myfile << "\n";
	}
	myfile.close();

	END_TRY_CATCH(;)
}