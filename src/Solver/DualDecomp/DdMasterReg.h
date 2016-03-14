/*
 * DdMasterReg.h
 *
 *  Created on: Feb 11, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DUALDECOMP_DDMASTERREG_H_
#define SRC_SOLVER_DUALDECOMP_DDMASTERREG_H_

#include "Solver/DualDecomp/DdMaster.h"

/**
 * Regularized bundle method
 */
class DdMasterReg: public DdMaster {
public:

	/** constructor */
	DdMasterReg(DspParams * par, DecModel * model, StoMessage * message, int nworkers, int maxnumsubprobs);

	/** destructor */
	virtual ~DdMasterReg();

	/** initialize */
	virtual STO_RTN_CODE init();

	/** solve */
	virtual STO_RTN_CODE solve();

	/** update problem */
	virtual STO_RTN_CODE updateProblem();

	/** set init solution */
	virtual STO_RTN_CODE setInitSolution(const double * sol);

protected:

	/** create problem */
	virtual STO_RTN_CODE createProblem();

private:

	/** update objective */
	STO_RTN_CODE refreshObjective();

	/** solve without regularization */
	STO_RTN_CODE terminationTest();

private:

	double sum_of_thetas_; /**< sum of thetas */
	int nthetas_;          /**< number of thetas */
	int nlambdas_;         /**< number of lambdas */

	double stability_param_;    /**< stability parameter */
	double * stability_center_; /**< stability center */

	int * nlastcuts_; /**< number of cuts generated at the last iteration */

	double ** primsol_to_worker_; /**< primal solution (theta and lambda) given to each worker */

	/** To reduce overhead for memory allocation */
	int * irowQ_;  /**< Hessian row index */
	int * jcolQ_;  /**< Hessian column index */
	double * dQ_;  /**< lower triangle Hessian elements */
	double * obj_; /**< linear objective coefficient */

	bool doSolve_; /**< indicate whether to solve master or to use the current solution stored. */
};

#endif /* SRC_SOLVER_DUALDECOMP_DDMASTERREG_H_ */
