/*
 * DdMasterAdmmLinearProximal.h
 *
 *  Created on: Apr 5, 2023
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DUALDECOMP_DDMASTERADMMLINEARPROXIMAL_H_
#define SRC_SOLVER_DUALDECOMP_DDMASTERADMMLINEARPROXIMAL_H_

#include "Solver/DualDecomp/DdMaster.h"

/** A class for implementing an inexact ADMM master solver using linear approximation. */
class DdMasterAdmmLinearProximal: public DdMaster {
public:

	/** A default constructor. */
	DdMasterAdmmLinearProximal(
			DecModel *   model,  /**< model pointer */
			DspParams *  par,    /**< parameter pointer */
			DspMessage * message /**< message pointer */);

	/** A copy constructor. */
	DdMasterAdmmLinearProximal(const DdMasterAdmmLinearProximal& rhs);

	/** A default destructor. */
	virtual ~DdMasterAdmmLinearProximal();

	/** A clone function */
	virtual DdMasterAdmmLinearProximal* clone() const {
		return new DdMasterAdmmLinearProximal(*this);
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
	double rho_;				/**< coefficient for quadratic term */
	double * gradient_;    		/**< subgradient */
	double * multipliers_; 		/**< Lagrangian multipliers */
	bool isPrimFeas;			/**< Feasibility of lagrangian multipliers */
    double * sum_mean_; 		/**< weighted sum of mean lagrangian multipliers */
	double * mean_multipliers_; /**< mean lagrangian multipliers */
	double auglagrange;			/**< augmented lagrangian function value */
	double primres;				/**< primal residual */
	double dualres;				/**< dual residual */

public:
	/** solver statistics */
	vector<double>		 s_auglagrange_;	/**< history of augmented lagrange function */
	vector<double>       s_primres_;  		/**< history of primal residual values */
	vector<double>		 s_dualres_;		/**< history of dual residual values*/
};

#endif /* SRC_SOLVER_DUALDECOMP_DDMASTERADMMLINEARPROXIMAL_H_ */
