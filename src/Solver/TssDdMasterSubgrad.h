/*
 * TssDdMasterSubgrad.h
 *
 *  Created on: Mar 3, 2015
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_TSSDDMASTERSUBGRAD_H_
#define SRC_SOLVER_TSSDDMASTERSUBGRAD_H_

#include "Solver/TssDdMaster.h"

class TssDdMasterSubgrad: public TssDdMaster
{
public:
	TssDdMasterSubgrad(StoParam * par) :
		TssDdMaster(par),
		ncols_(0),
		ncols_first_(0),
		nscenarios_(0),
		nstalls_(0),
		stepscal_(2.),
		stepsize_(0.),
		gradient_(NULL),
		multipliers_(NULL) {}

	virtual ~TssDdMasterSubgrad();

	/** create problem */
	virtual STO_RTN_CODE createProblem(const TssModel * model);

	/** update problem: may update dual bound */
	virtual STO_RTN_CODE updateProblem(
			double primal_bound, /**< primal bound of the original problem */
			double & dual_bound, /**< dual bound of the original problem */
			double * objvals,    /**< objective values of subproblems */
			double ** solution   /**< subproblem solutions */);

	/** solve problem */
	virtual STO_RTN_CODE solve();

	/** solution status */
	virtual STO_RTN_CODE getStatus() {return STO_STAT_OPTIMAL;}

	/** get Lagrangian multiplier */
	virtual const double * getLagrangian() {return solution_;}

	double getConstScalar() {return stepscal_;}

private:

	int ncols_;            /**< dimension of the subgradient */
	int ncols_first_;      /**< number of first-stage columns */
	int nscenarios_;       /**< number of scenarios */
	int nstalls_;          /**< number of iterations UB were not improved */
	double stepscal_;      /**< another scalar between 0 and 2 */
	double stepsize_;      /**< stepsize */
	double * gradient_;    /**< subgradient */
	double * multipliers_; /**< Lagrangian multipliers */
};

#endif /* SRC_SOLVER_TSSDDMASTERSUBGRAD_H_ */
