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

	/** write output to a file */
	virtual void write(const char * filename);

protected:

	/** create problem */
	virtual DSP_RTN_CODE createProblem();

private:

	int itercount_;		   		/**< number of iterations */
	int nstalls_;          		/**< number of iterations UB were not improved */
	double stepscal_;      		/**< another scalar between 0 and 2 */
	double stepsize_;      		/**< stepsize */
	double * gradient_;    		/**< subgradient */
	double * multipliers_; 		/**< Lagrangian multipliers */
	bool isPrimFeas;			/**< Feasibility of lagrangian multipliers */
    double * prev_gradient_; 	/**< mean subgradient at previous iteration */
	double * mean_multipliers_; /**< mean lagrangian multipliers */

public:
	/** solver statistics */
	vector<double>       s_primres_;  /**< history of primal residual values */

};

#endif /* SRC_SOLVER_DUALDECOMP_DDMASTERADMMLINEAR_H_ */
