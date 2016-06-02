/*
 * DecDdMasterSubgrad.h
 *
 *  Created on: Mar 3, 2015
 *      Author: kibaekkim, ctjandra
 */

#ifndef SRC_SOLVER_DECDDMASTERSUBGRAD_H_
#define SRC_SOLVER_DECDDMASTERSUBGRAD_H_

#include "Solver/DecDdMaster.h"

class DecDdMasterSubgrad: public DecDdMaster
{
public:
	DecDdMasterSubgrad(DspParams * par) :
		DecDdMaster(par),
		ncoupling_(0),
		nsubprobs_(0),
		nstalls_(0),
		stepscal_(2.),
		stepsize_(0.),
		gradient_(NULL),
		multipliers_(NULL),
		nonanticipativity_(false),
		model_(NULL) {}

	virtual ~DecDdMasterSubgrad();

	/** create problem */
	virtual DSP_RTN_CODE createProblem(DecModel * model);

	/** update problem: may update dual bound */
	virtual DSP_RTN_CODE updateProblem(
			double primal_bound,     /**< primal bound of the original problem */
			double & dual_bound,     /**< dual bound of the original problem */
			double * primal_objvals, /**< objective values of subproblems */
			double * dual_objvals,   /**< objective values of subproblems */
			double ** solution       /**< subproblem solutions */);

	/** solve problem */
	virtual DSP_RTN_CODE solve();

	/** solution status */
	virtual DSP_RTN_CODE getStatus() {return STO_STAT_OPTIMAL;}

	/** get Lagrangian multiplier */
	virtual const double * getLagrangian() {return solution_;}

	/** termination test */
	virtual bool terminate();

	double getConstScalar() {return stepscal_;}

private:

	int ncoupling_;        /**< number of coupling constraints (i.e., dimension of Lagrange multipliers and subgradient) */
	int nsubprobs_;        /**< number of subproblems */
	int nstalls_;          /**< number of iterations UB were not improved */
	double stepscal_;      /**< another scalar between 0 and 2 */
	double stepsize_;      /**< stepsize */
	double * gradient_;    /**< subgradient */
	double * multipliers_; /**< Lagrangian multipliers */
	bool nonanticipativity_; /**< whether the constraints are nonanticipativity constraints */
	DecModel * model_;     /**< model */
};

#endif /* SRC_SOLVER_DECDDMASTERSUBGRAD_H_ */
