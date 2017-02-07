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
#include "SimpleVector.h"

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
		const double gamma = 1.0;
		bool wellcentered = true;
		double objval     = qpdata->objectiveValue(qpvars);
		double dualityGap = resids->dualityGap();
		double relDualityGap = dualityGap / (1. + fabs(objval));

		double* z = dynamic_cast<SimpleVector*>(qpvars->z.ptr())->elements();
		double* s = dynamic_cast<SimpleVector*>(qpvars->s.ptr())->elements();
		int mz = qpdata->mz;
		double sz;
		for (int i = 0; i < mz; ++i) {
			sz = s[i] * z[i];
			if (sz < mu * gamma || sz > mu / gamma) {
				wellcentered = false;
				break;
			}
		}
#ifdef DSP_DEBUG
		printf("OoqpStatus: mu %e, gamma %e, objval %e, lowerBound %e, dualityGap %e (%e), epsilon %e\n",
				mu, gamma, objval, lowerBound_, dualityGap, relDualityGap, epsilon_);
#endif

		if (wellcentered > 0) {
			if (objval <= lowerBound_)
				return 0; /**< successful termination */
			if (relDualityGap < epsilon_)
				return 0; /**< successful termination */
#if 0
			/** What did I try to do here? */
			if (relDualityGap < epsilon_ &&
					upperBound_ < COIN_DBL_MAX &&
					upperBound_ - objval > 1.e-4 * fabs(upperBound_))
				return 0; /**< successful termination */
#endif
		}

#if 0
		/** What did I try to do here? */
		if (mu <= solver->getMuTol() &&
				resids->residualNorm() <= solver->getArTol() * data->datanorm())
			return 0; /**< successful termination */
#endif
	}

	return solver->defaultStatus(data, vars, resids, i, mu, level);
}
