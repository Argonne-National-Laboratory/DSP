/*
 * BdDriver.h
 *
 *  Created on: Feb 17, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_BENDERS_BDDRIVER_H_
#define SRC_SOLVER_BENDERS_BDDRIVER_H_

#include "Solver/DspDriver.h"
#include "Solver/Benders/BdMW.h"

class BdDriver: public DspDriver {

public:

	/** constructor */
	BdDriver(
			DspParams * par, /**< parameter pointer */
			DecModel * model /**< model pointer */);

	/** destructor */
	virtual ~BdDriver();

	/** initialize */
	virtual DSP_RTN_CODE init() {return DSP_RTN_OK;}

	/** run */
	virtual DSP_RTN_CODE run() {return DSP_RTN_OK;}

	/** finalize */
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
