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
#include "SolverInterface/SolverInterfaceOoqp.h"

class OoqpEps: public SolverInterfaceOoqp
{
public:

	/** default constructor */
	OoqpEps(DspParams * par) :
		SolverInterfaceOoqp(par),
		hasOoqpStatus_(true),
		epsilon_(COIN_DBL_MAX),
		lowerBound_(-COIN_DBL_MAX),
		upperBound_(COIN_DBL_MAX),
		dualityGap_(0.0),
		suboptimal_(false) {}

	/** default destructor */
	virtual ~OoqpEps() {}

	/** solve */
	virtual void solve();

	/** set my OOQP status */
	void setOoqpStatus(double epsilon, double lowerBound, double upperBound);

	/** core part for load problem */
//	virtual void gutsOfLoadProblem();

	/** get duality gap tolerance */
	double getDualityGapTol() {return epsilon_;}

	/** get relative duality gap */
	double getDualityGap() {return dualityGap_;}

	/** is suboptimal? */
	bool isSuboptimal() {return suboptimal_;}

public:

	bool hasOoqpStatus_; /**< has user-defined OOQP status object? */

protected:

	double epsilon_;
	double lowerBound_;
	double upperBound_;
	double dualityGap_; /** duality gap */
	bool suboptimal_; /**< indicating whether the solution is suboptimal */
};

#endif /* SRC_SOLVER_OOQPEPS_H_ */
