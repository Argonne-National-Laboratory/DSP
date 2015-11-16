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
#include "SolverInterface/SolverInterface.h"
#include "Utility/DspParams.h"
#include "Model/DecModel.h"

using namespace std;

class DecDdMaster
{
public:

	/** default constructor */
	DecDdMaster(DspParams * par) :
		par_(par),
		objval_(COIN_DBL_MAX),
		solution_(NULL),
		parLogLevel_(par->getIntParam("LOG_LEVEL")),
		parTR_(par->getBoolParam("DD/TR")),
		parTrDecrease_(par->getBoolParam("DD/TR/DECREASE")),
		parMasterAlgo_(par->getIntParam("DD/MASTER_ALGO")),
		parStopTol_(par->getDblParam("DD/STOP_TOL")),
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
			double primal_bound,     /**< primal bound of the original problem */
			double & dual_bound,     /**< dual bound of the original problem */
			double * primal_objvals, /**< objective values of subproblems */
			double * dual_objvals,   /**< objective values of subproblems */
			double ** solution       /**< subproblem solutions */) = 0;

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

	/** termination test */
	virtual bool terminate() = 0;

protected:

	DspParams * par_;   /**< parameters */

	double objval_;
	double * solution_; /**< solution */

	/** parameters */
	int parLogLevel_;
	bool parTR_;
	bool parTrDecrease_;
	int parMasterAlgo_;
	double parStopTol_;

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
