/*
 * OoqpEps.cpp
 *
 *  Created on: Jan 28, 2015
 *      Author: kibaekkim
 */

/** DSP */
#include <Utility/DspMacros.h>
#include "SolverInterface/OoqpEps.h"
#include "SolverInterface/OoqpStatus.h"

/** set my OOQP status */
void OoqpEps::setOoqpStatus(double epsilon, double lowerBound, double upperBound)
{
	releaseOOQP();
	epsilon_       = epsilon;
	lowerBound_    = lowerBound;
	upperBound_    = upperBound;
}

/** solve */
void OoqpEps::solve()
{
	if (released_)
	{
		GUTS_OF_LOAD_PROBLEM();
		if (hasOoqpStatus_)
		{
			OoqpStatus * mystat = new OoqpStatus(epsilon_, lowerBound_, upperBound_);
			solver_->useStatus(mystat);
		}
	}

	//prob_->print();
	if (print_level_ > 0)
		solver_->monitorSelf();
	int status = solver_->solve(prob_, vars_, resid_);
	switch(status)
	{
	case 0:
	{
		status_ = DSP_STAT_OPTIMAL;
		objval_ = prob_->objectiveValue(vars_);

		/** get variable values */
		vars_->x->copyIntoArray(x_);
		vars_->y->copyIntoArray(y_);
		vars_->lambda->copyIntoArray(lambda_);
		vars_->pi->copyIntoArray(pi_);
		vars_->gamma->copyIntoArray(gamma_);
		vars_->phi->copyIntoArray(phi_);

		nLambdas_ =  vars_->lambda->length();
		nPis_ = vars_->pi->length();
		dualityGap_ = resid_->dualityGap();

		/** number of iterations */
		nIters_ = solver_->iter;

		/** suboptimal? */
		if (vars_->mu() <= solver_->getMuTol() &&
			resid_->residualNorm() <= solver_->getArTol() * prob_->datanorm())
			suboptimal_ = false;
		else
			suboptimal_ = true;
		//printf("optimal? %s\n", suboptimal_ ? "no" : "yes");

		break;
	}
	case 1:
		status_ = DSP_STAT_STOPPED_UNKNOWN;
		break;
	case 2:
		status_ = DSP_STAT_STOPPED_ITER;
		break;
	case 3:
		status_ = DSP_STAT_PRIM_INFEASIBLE;
		break;
	case 4:
	default:
		status_ = DSP_STAT_UNKNOWN;
		break;
	}
}
