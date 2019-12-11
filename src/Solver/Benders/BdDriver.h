/*
 * BdDriver.h
 *
 *  Created on: Feb 17, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_BENDERS_BDDRIVER_H_
#define SRC_SOLVER_BENDERS_BDDRIVER_H_

#include "Solver/DecSolver.h"
#include "Solver/Benders/BdMW.h"

/** A base driver class for Benders decomposition */
class BdDriver: public DecSolver {

public:

	/** A default constructor. */
	BdDriver(
			DecModel *   model,  /**< model pointer */
			DspParams *  par,    /**< parameters */
			DspMessage * message /**< message pointer */);

	/** A copy constructor. */
	BdDriver(const BdDriver&rhs);

	/** A default destructor. */
	virtual ~BdDriver();

	/** A clone function */
	virtual BdDriver* clone() const {
		return new BdDriver(*this);
	}

	/** A virtual member for initializing the driver. */
	virtual DSP_RTN_CODE init() {return DSP_RTN_OK;}

	/** A virtual member for running the driver. */
	virtual DSP_RTN_CODE run() {return solve();}
	virtual DSP_RTN_CODE solve() {return DSP_RTN_ERR;}

	/** A virtual memeber for finalizing the driver. */
	virtual DSP_RTN_CODE finalize() {return DSP_RTN_OK;}

protected:

	/** find lower bound */
	virtual DSP_RTN_CODE findLowerBound() {return DSP_RTN_OK;}

	/** collect second-stage solutions */
	virtual DSP_RTN_CODE collectSolution() {return DSP_RTN_OK;}

public:

	/** set dual objective */
	virtual void setDualObjective(double dualobj) {dualobj_ = dualobj;}

	/** set auxiliary variable data */
	virtual DSP_RTN_CODE setAuxVarData(
			int      size, /**< size of arrays */
			double * obj,  /**< objective function coefficients */
			double * clbd, /**< column lower bounds */
			double * cubd  /**< column upper bounds */);

	/** set initial solution */
	virtual DSP_RTN_CODE setSolution(
			int      size,    /**< size of array */
			double * solution /**< solution */);

	/** set branch priorities */
	virtual DSP_RTN_CODE setPriorities(
			int   size,      /**< size of array */
			int * priorities /**< branch priority */);

protected:

	BdMW * mw_; /**< master-worker manager */

	int      aux_size_; /**< number of auxiliary variables */
	double * aux_obj_;  /**< objective coefficients */
	double * aux_clbd_; /**< column lower bound */
	double * aux_cubd_; /**< column upper bound */

	Solutions initsols_; /**< set of initial solutions */
	int numPriorities_;  /**< length of branch priorities */
	int * priorities_;   /**< branch priority */
};

#endif /* SRC_SOLVER_BENDERS_BDDRIVER_H_ */
