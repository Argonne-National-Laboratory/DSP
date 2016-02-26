/*
 * DspDriver.h
 *
 *  Created on: Feb 16, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DSPDRIVER_H_
#define SRC_SOLVER_DSPDRIVER_H_

#include "Utility/DspParams.h"
#include "Utility/StoMessage.h"
#include "Model/DecModel.h"

class DspDriver {
public:

	/** constructor */
	DspDriver(DspParams * par, DecModel * model);

	/** destructor */
	virtual ~DspDriver();

	/** initialize */
	virtual STO_RTN_CODE init() = 0;

	/** run */
	virtual STO_RTN_CODE run() = 0;

public:

	/** get solution status */
	virtual STO_RTN_CODE getStatus() {return status_;}

	/** get primal solution */
	virtual const double * getPrimalSolution() {return primsol_;}

	/** get dual solution */
	virtual const double * getDualSolution() {return dualsol_;}

	/** get primal objective value */
	virtual double getPrimalObjectiveValue() {return primobj_;}

	/** get dual objective value */
	virtual double getDualObjectiveValue() {return dualobj_;}

public:

	/** get solution time */
	virtual double getCpuTime() {return cputime_;}
	virtual double getWallTime() {return walltime_;}

	/** get number of iterations */
	virtual int getNumIterations() {return numIterations_;}

	/** get number of branch-and-bound nodes */
	virtual int getNumNodes() {return numNodes_;}

protected:

	DspParams * par_;
	DecModel * model_;
	StoMessage * message_; /**< message handler */

protected:

	STO_RTN_CODE status_; /**< solution status */
	double * primsol_;    /**< primal solution in extensive form */
	double * dualsol_;    /**< dual solution in extensive form */
	double   primobj_;    /**< primal objective bound */
	double   dualobj_;    /**< dual objective bound */

protected:

	double cputime_;    /**< cpu time */
	double walltime_;   /**< wall time */
	int numIterations_; /**< number of iterations for a given method */
	int numNodes_;      /**< number of branch-and-bound tree nodes */
};

#endif /* SRC_SOLVER_DSPDRIVER_H_ */
