/*
 * TssSolver.h
 *
 *  Created on: Sep 24, 2014
 *      Author: kibaekkim
 */

#ifndef TSSSOLVER_H_
#define TSSSOLVER_H_

/** DSP */
#include "Utility/StoConfig.h"
#include "Utility/StoRtnCodes.h"
#include "Utility/StoUtility.h"
#include "Utility/StoMessage.h"
#include "Utility/DspParams.h"
#include "Model/TssModel.h"

class TssSolver {
public:

	/** default constructor */
	TssSolver();

	/** default destructor */
	virtual ~TssSolver();

	/** load model object; will have a shallow pointer (not deep copy) */
	STO_RTN_CODE loadModel(DspParams * par, TssModel * model);

	/** solve */
	virtual STO_RTN_CODE solve() = 0;

	/** get message handler */
	CoinMessageHandler * messageHandler() {return handler_;}

protected:

	/**
	 * This converts ctype to the number of integer variables and their indices.
	 */
	void convertColTypes(int ncols, char * ctype, int & len, int *& ind);

protected:

	TssModel * model_; /**< TssModel object */
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

#endif /* TSSSOLVER_H_ */
