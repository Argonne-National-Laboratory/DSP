/*
 * DecSolver.h
 *
 *  Created on: Sep 24, 2014
 *      Author: kibaekkim, ctjandra
 */

#ifndef DECSOLVER_H_
#define DECSOLVER_H_

/** DSP */
#include <Utility/DspMessage.h>
#include <Utility/DspMpi.h>
#include <Utility/DspRtnCodes.h>
#include "Utility/StoConfig.h"
#include "Utility/DspParams.h"
#include "Model/DecModel.h"

class DecSolver {
public:

	/** default constructor */
	DecSolver();

	/** default destructor */
	virtual ~DecSolver();

	/** load model object; will have a shallow pointer (not deep copy) */
	virtual DSP_RTN_CODE loadModel(DspParams * par, DecModel * model);

	/** solve */
	virtual DSP_RTN_CODE solve() = 0;

protected:

	DecModel * model_;     /**< DecModel object */
	DspParams * par_;      /**< parameters */
	DspMessage * message_; /**< message */

	double time_remains_; /**< wall clock time remains */

	/** parameters */
	int parLogLevel_;
	int parNodeLim_;
	int parIterLim_;
	double parWallLim_;
	int parScipDispFreq_;
	double parScipGapTol_;
	double parScipTimeLim_;
	const bool * parRelaxIntegrality_;

public:

	/** global clocks */
	double ctime_start_; /**< start time (cpu) */
	double wtime_start_; /**< start time (wall) */

	/**
	 * Solutions and Statistics
	 */

	DSP_RTN_CODE status_; /**< solution status */
	double * solution_;   /**< solution in extensive form */
	double primalBound_;  /**< primal objective bound */
	double dualBound_;    /**< dual objective bound */

	int numIterations_;   /**< number of iterations for a given method */
	int numNodes_;        /**< number of branch-and-bound tree nodes */
	int clockType_;       /**< 0: cpu, 1: wall */
	double solutionTime_; /**< solution time */
};

#endif /* DECSOLVER_H_ */
