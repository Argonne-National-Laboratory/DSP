/*
 * DecSolver.h
 *
 *  Created on: Sep 24, 2014
 *      Author: kibaekkim, ctjandra
 */

#ifndef DECSOLVER_H_
#define DECSOLVER_H_

/** DSP */
#include "Utility/StoConfig.h"
#include "Utility/StoRtnCodes.h"
#include "Utility/StoUtility.h"
#include "Utility/StoMessage.h"
#include "Utility/DspParams.h"
#include "Model/DecModel.h"
//#include "Solver/StoParam.h"

class DecSolver {
public:

	/** default constructor */
	DecSolver();

	/** default destructor */
	virtual ~DecSolver();

	/** load model object; will have a shallow pointer (not deep copy) */
	virtual STO_RTN_CODE loadModel(DspParams * par, DecModel * model);

	/** solve */
	virtual STO_RTN_CODE solve() = 0;

	/** get message handler */
	CoinMessageHandler * messageHandler() {return handler_;}

protected:

	DecModel * model_; /**< DecModel object */
	DspParams * par_;   /**< parameters */

	CoinMessageHandler * handler_; /**< message handler */
	StoMessage * message_;         /**< message */

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

	/**
	 * Solutions and Statistics
	 */

	STO_RTN_CODE status_; /**< solution status */
	double * solution_;   /**< solution in extensive form */
	double primalBound_;  /**< primal objective bound */
	double dualBound_;    /**< dual objective bound */

	int numIterations_;   /**< number of iterations for a given method */
	int numNodes_;        /**< number of branch-and-bound tree nodes */
	int clockType_;       /**< 0: cpu, 1: wall */
	double solutionTime_; /**< solution time */
};

#endif /* DECSOLVER_H_ */
