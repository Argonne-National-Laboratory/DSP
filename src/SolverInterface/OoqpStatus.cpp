/*
 * OoqpStatus.cpp
 *
 *  Created on: Jan 23, 2015
 *      Author: kibaekkim
 */

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
#include "Utility/DspMessage.h"
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
	if (qpdata && qpvars && epsilon_ > 0)
	{
		const double gamma = 0.5;
		bool wellcentered = true;

		/** retrieve OOQP data and variables */
		int mz = qpdata->mz;
		double* iclow = dynamic_cast<SimpleVector*>(qpdata->iclow.ptr())->elements();
		double* icupp = dynamic_cast<SimpleVector*>(qpdata->icupp.ptr())->elements();
		double* t = dynamic_cast<SimpleVector*>(qpvars->t.ptr())->elements();
		double* u = dynamic_cast<SimpleVector*>(qpvars->u.ptr())->elements();
		double* lambda = dynamic_cast<SimpleVector*>(qpvars->lambda.ptr())->elements();
		double* pi = dynamic_cast<SimpleVector*>(qpvars->pi.ptr())->elements();

		double sz;
		for (int i = 0; i < mz; ++i) {
			sz = 0.0;
			if (iclow[i] > 0)
				sz += t[i]*lambda[i];
			if (icupp[i] > 0)
				sz += u[i]*pi[i];
			if (sz < mu * gamma || sz > mu / gamma) {
				wellcentered = false;
				break;
			}
		}

#ifdef DSP_DEBUG_MORE
		printf("OoqpStatus: mu %e, gamma %e, objval %e, lowerBound %e, dualityGap %e (%e), epsilon %e, %s\n",
				mu, gamma, objval, lowerBound_, dualityGap, relDualityGap, epsilon_, wellcentered ? "centered" : "not");
#endif

		if (wellcentered) {
			double objval     = qpdata->objectiveValue(qpvars);
			double dualityGap = resids->dualityGap();
			double relDualityGap = dualityGap / (1. + fabs(objval));
//			if (objval <= lowerBound_)
//				return 0; /**< successful termination */
			if (dualityGap > 0 && relDualityGap < epsilon_) {
				DSPdebugMessage("relative duality gap %e < epsilon %e\n", relDualityGap, epsilon_);
				return 0; /**< successful termination */
			}
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
