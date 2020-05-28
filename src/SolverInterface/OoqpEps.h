/*
 * OoqpEps.h
 *
 *  Created on: Jan 28, 2015
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_OOQPEPS_H_
#define SRC_SOLVER_OOQPEPS_H_

/** DSP */
#include "SolverInterface/OoqpStatus.h"
#include "SolverInterface/OsiOoqpSolverInterface.hpp"

class OoqpEps: public OsiOoqpSolverInterface
{
public:

	/** default constructor */
	OoqpEps() :
		OsiOoqpSolverInterface(),
		hasOoqpStatus_(true),
		mystat_(NULL),
		epsilon_(COIN_DBL_MAX),
		lowerBound_(-COIN_DBL_MAX),
		upperBound_(COIN_DBL_MAX),
		dualityGap_(0.0),
		suboptimal_(false) {}

	/** default destructor */
	virtual ~OoqpEps() {
		delete mystat_;
		mystat_ = NULL;
	}

	/** solve */
	virtual void resolve();

	/** set my OOQP status */
	void setOoqpStatus(double epsilon, double lowerBound, double upperBound);

	/** get duality gap tolerance */
	double getDualityGapTol() {return epsilon_;}

	/** get relative duality gap */
	double getDualityGap() {return dualityGap_;}

	/** is suboptimal? */
	bool isSuboptimal() {return suboptimal_;}

public:

	bool hasOoqpStatus_; /**< has user-defined OOQP status object? */

protected:

	OoqpStatus* mystat_; /**< my customized status for early termination */
	double epsilon_; /**< duality gap tolerance */
	double lowerBound_;
	double upperBound_;
	double dualityGap_; /** duality gap */
	bool suboptimal_; /**< indicating whether the solution is suboptimal */
};

#endif /* SRC_SOLVER_OOQPEPS_H_ */
