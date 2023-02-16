/*
 * DdMasterAdmmLinear.h
 *
 *  Created on: Feb 10, 2023
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DUALDECOMP_DDMASTERADMMLINEAR_H_
#define SRC_SOLVER_DUALDECOMP_DDMASTERADMMLINEAR_H_

#include "Solver/DualDecomp/DdMaster.h"

/** A class for implementing an inexact ADMM master solver using linear approximation. */
class DdMasterAdmmLinear: public DdMaster {
public:

	/** A default constructor. */
	DdMasterAdmmLinear(
			DecModel *   model,  /**< model pointer */
			DspParams *  par,    /**< parameter pointer */
			DspMessage * message /**< message pointer */);

	/** A copy constructor. */
	DdMasterAdmmLinear(const DdMasterAdmmLinear& rhs);

	/** A default destructor. */
	virtual ~DdMasterAdmmLinear();

	/** A clone function */
	virtual DdMasterAdmmLinear* clone() const {
		return new DdMasterAdmmLinear(*this);
	}

	/** A virtual member for initializing solver. */
	virtual DSP_RTN_CODE init();

	/** A virtual member for solving problem. */
	virtual DSP_RTN_CODE solve();

	/** update problem */
	virtual DSP_RTN_CODE updateProblem();

	/** get number of lambdas */
	virtual int getNumLambdas() {return model_->getNumCouplingRows();}

protected:

	/** create problem */
	virtual DSP_RTN_CODE createProblem();

private:

	int nstalls_;          /**< number of iterations UB were not improved */
	double stepscal_;      /**< another scalar between 0 and 2 */
	double stepsize_;      /**< stepsize */
	double * gradient_;    /**< subgradient */
    double * prev_gradient_; /** subgradient at previous iteration */
	double * multipliers_; /**< Lagrangian multipliers */
};

#endif /* SRC_SOLVER_DUALDECOMP_DDMASTERADMMLINEAR_H_ */
