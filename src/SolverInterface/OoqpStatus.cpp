/*
 * OoqpStatus.cpp
 *
 *  Created on: Jan 23, 2015
 *      Author: kibaekkim
 */

#include <Utility/DspMessage.h>
#include "stdio.h"
#include "math.h"

//#define DSP_DEBUG

/** OOQP */
#include "QpGenData.h"
#include "QpGenVars.h"
#include "QpGenResiduals.h"
#include "GondzioSolver.h"
#include "MehrotraSolver.h"

/** DSP */
#include "SolverInterface/OoqpStatus.h"

OoqpStatus::~OoqpStatus()
{
	// TODO Auto-generated destructor stub
}

int OoqpStatus::doIt(
		Solver *    solver,
		Data *      data,
		Variables * vars,
		Residuals * resids,
		int         i,
		double      mu,
		int         level)
{
	QpGenData * qpdata = dynamic_cast<QpGenData*>(data);
	QpGenVars * qpvars = dynamic_cast<QpGenVars*>(vars);
	if (qpdata && qpvars)
	{
		double objval     = qpdata->objectiveValue(qpvars);
		double dualityGap = resids->dualityGap();
		double relDualityGap = dualityGap / (1. + fabs(objval));
		double mc = dualityGap / mu;
#ifdef DSP_DEBUG
		printf("KK: mu %e rnorm %e dnorm %e Ar %e dualityGap %e objval %e",
				mu, resids->residualNorm(), data->datanorm(),
				resids->residualNorm() / data->datanorm(), dualityGap, objval);
		if (lowerBound_ > -COIN_DBL_MAX)
			printf(" lowerBound %e", lowerBound_);
		printf(" upperBound %e epsilon %e\n",
				upperBound_, epsilon_);
#endif
		if (mc > 0)
		{
			if (objval <= lowerBound_)
				return 0; /**< successful termination */
			if (relDualityGap < epsilon_ &&
					upperBound_ < COIN_DBL_MAX &&
					upperBound_ - objval > 1.e-4 * fabs(upperBound_))
				return 0; /**< successful termination */
		}
		if (mu <= solver->getMuTol() &&
				resids->residualNorm() <= solver->getArTol() * data->datanorm())
			return 0; /**< successful termination */
	}

	return solver->defaultStatus(data, vars, resids, i, mu, level);
}
