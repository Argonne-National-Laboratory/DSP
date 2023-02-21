/*
 * DdMasterSubgrad.h
 *
 *  Created on: Feb 16, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DUALDECOMP_DDMASTERSUBGRAD_H_
#define SRC_SOLVER_DUALDECOMP_DDMASTERSUBGRAD_H_

#include "Solver/DualDecomp/DdMaster.h"

/** A class for implementing a subgradient master solver. */
class DdMasterSubgrad: public DdMaster {
public:

	/** A default constructor. */
	DdMasterSubgrad(
			DecModel *   model,  /**< model pointer */
			DspParams *  par,    /**< parameter pointer */
			DspMessage * message /**< message pointer */);

	/** A copy constructor. */
	DdMasterSubgrad(const DdMasterSubgrad& rhs);

	/** A default destructor. */
	virtual ~DdMasterSubgrad();

	/** A clone function */
	virtual DdMasterSubgrad* clone() const {
		return new DdMasterSubgrad(*this);
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

	int itercount_;		   /**< number of iterations */
	int nstalls_;          /**< number of iterations UB were not improved */
	double stepscal_;      /**< another scalar between 0 and 2 */
	double stepsize_;      /**< stepsize */
	double * gradient_;    /**< subgradient */
	double * multipliers_; /**< Lagrangian multipliers */
};

#endif /* SRC_SOLVER_DUALDECOMP_DDMASTERSUBGRAD_H_ */
