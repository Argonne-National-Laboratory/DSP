/*
 * OoqpEps.cpp
 *
 *  Created on: Jan 28, 2015
 *      Author: kibaekkim
 */

//#define DSP_DEBUG

/** DSP */
#include "Utility/DspMacros.h"
#include "Utility/DspMessage.h"
#include "SolverInterface/OoqpEps.h"
#include "SolverInterface/OoqpStatus.h"
/** Ooqp */
#include "GondzioSolver.h"
#include "MehrotraSolver.h"

/** set my OOQP status */
void OoqpEps::setOoqpStatus(double epsilon, double lowerBound, double upperBound) {
	epsilon_       = epsilon;
	lowerBound_    = lowerBound;
	upperBound_    = upperBound;
}

/** solve */
void OoqpEps::resolve() {
	if (updated_) {
		/** free Ooqp objects */
		freeCachedOoqp();

		/** convert Osi data to Ooqp data */
		convertOsiToOoqp(qpgen_, prob_);
		DSPdebug(prob_->print());

		/** declare variables */
		vars_ = (QpGenVars*)qpgen_->makeVariables(prob_);

		/** declare residuals */
		resid_ = (QpGenResiduals*)qpgen_->makeResiduals(prob_);

		/** create solver */
		//solver_ = new GondzioSolver(qp_, prob_);
		solver_ = new MehrotraSolver(qpgen_, prob_);

		/** assign my status */
		if (hasOoqpStatus_) {
			if (mystat_) {
				delete mystat_;
				mystat_ = NULL;
			}
			mystat_ = new OoqpStatus(epsilon_, lowerBound_, upperBound_);
			solver_->useStatus(mystat_);
		}
	} else {
		mystat_->setStatus(epsilon_, lowerBound_, upperBound_);
	}

	/** actual solution */
	gutsOfSolve();

	/** duality gap */
	dualityGap_ = resid_->dualityGap();

	/** is suboptimal? */
	if (vars_->mu() <= solver_->getMuTol() &&
			resid_->residualNorm() <= solver_->getArTol() * prob_->datanorm()) {
		suboptimal_ = false;
	} else
		suboptimal_ = true;
}
