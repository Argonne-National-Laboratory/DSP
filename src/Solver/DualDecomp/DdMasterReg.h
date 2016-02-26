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

protected:

	/** create problem */
	virtual STO_RTN_CODE createProblem();

private:

	/** update objective */
	STO_RTN_CODE refreshObjective();

	/** solve without regularization */
	STO_RTN_CODE solveWithoutReg();

private:

	int nthetas_;  /**< number of thetas */
	int nlambdas_; /**< number of lambdas */

	double stability_param_;    /**< stability parameter */
	double * stability_center_; /**< stability center */

	int * nlastcuts_; /**< number of cuts generated at the last iteration */

	/** To reduce overhead for memory allocation */
	int * irowQ_;  /**< Hessian row index */
	int * jcolQ_;  /**< Hessian column index */
	double * dQ_;  /**< lower triangle Hessian elements */
	double * obj_; /**< linear objective coefficient */
};

#endif /* SRC_SOLVER_DUALDECOMP_DDMASTERREG_H_ */
