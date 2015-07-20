/*
 * OoqpStatus.h
 *
 *  Created on: Jan 23, 2015
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_OOQPSTATUS_H_
#define SRC_SOLVER_OOQPSTATUS_H_

/** COIN */
#include "CoinFinite.hpp"

/** OOQP */
#include "Status.h"

class OoqpStatus: public Status
{
public:

	/** default constructor */
	OoqpStatus() :
		epsilon_(COIN_DBL_MAX),
		lowerBound_(-COIN_DBL_MAX),
		upperBound_(COIN_DBL_MAX) {}

	/** constructor with arguments */
	OoqpStatus(double epsilon, double lowerBound, double upperBound) :
		epsilon_(epsilon),
		lowerBound_(lowerBound),
		upperBound_(upperBound) {}

	virtual ~OoqpStatus();

	virtual int doIt(
			Solver *    solver,
			Data *      data,
			Variables * vars,
			Residuals * resids,
			int         i,
			double      mu,
			int         level);

private:

	double epsilon_;
	double lowerBound_;
	double upperBound_;
};

#endif /* SRC_SOLVER_OOQPSTATUS_H_ */
