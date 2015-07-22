/*
 * DecDdMaster.h
 *
 *  Created on: Dec 10, 2014
 *      Author: kibaekkim, ctjandra
 */

#ifndef SRC_SOLVER_DECDDMASTER_H_
#define SRC_SOLVER_DECDDMASTER_H_

#include <vector>
#include <iostream>

/** Coin-OR */
#include "OsiSolverInterface.hpp"
#include "OsiCuts.hpp"

/** DSP */
#include "Solver/SolverInterface.h"
#include "Solver/StoParam.h"
#include "Model/DecModel.h"

using namespace std;

class DecDdMaster
{
public:

	/** default constructor */
	DecDdMaster(StoParam * par) :
		par_(par),
		objval_(COIN_DBL_MAX),
		solution_(NULL),
		itncnt_accum_(0),
		wtime_accum_(0.)
	{
		/** nothing to do */
	}

	/** default destructor */
	virtual ~DecDdMaster() {}

	/** create problem */
	virtual STO_RTN_CODE createProblem(DecModel * model) = 0;

	/** update problem: may update dual bound */
	virtual STO_RTN_CODE updateProblem(
			double primal_bound, /**< primal bound of the original problem */
			double & dual_bound, /**< dual bound of the original problem */
			double * objvals,    /**< objective values of subproblems */
			double ** solution   /**< subproblem solutions */) = 0;

	/** solve problem */
	virtual STO_RTN_CODE solve() = 0;

	/** solution status */
	virtual STO_RTN_CODE getStatus() = 0;

	/** get Lagrangian multiplier */
	virtual const double * getLagrangian() = 0;

	/** get solution */
	virtual const double * getSolution() {return solution_;}

	/** get objective value */
	virtual double getObjValue() {return objval_;}

protected:

	StoParam * par_;   /**< parameters */

	double objval_;
	double * solution_; /**< solution */

public:

	vector<int>    itncnt_;         /**< total number of iterations */
	int            itncnt_accum_;

	vector<double> objvals_;
	vector<double> dual_objvals_;
	vector<double> primal_bounds_;
	vector<double> dual_bounds_;
	vector<double> wtime_solve_;
	double wtime_accum_;
};

#endif /* SRC_SOLVER_DECDDMASTER_H_ */
