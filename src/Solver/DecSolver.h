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

class DecSolver {
public:

	/** constructor */
	DecSolver(DspParams * par, DecModel * model, StoMessage * message);

	/** default destructor */
	virtual ~DecSolver();

	/** initialize */
	virtual STO_RTN_CODE init() = 0;

	/** finalize */
	virtual STO_RTN_CODE finalize() {return STO_RTN_OK;}

	/** solve */
	virtual STO_RTN_CODE solve() = 0;

public:

	/** solver status */
	virtual STO_RTN_CODE getStatus() {return status_;}

	/** get primal solution */
	virtual const double * getPrimalSolution() {return primsol_;}

	/** get dual solution */
	virtual const double * getDualSolution() {return dualsol_;}

	/** get primal objective */
	virtual double getPrimalObjective() {return primobj_;}

	/** get dual objective */
	virtual double getDualObjective() {return dualobj_;}

	/** get model pointer */
	virtual DecModel * getModelPtr() {return model_;}

	/** get parameter pointer */
	virtual DspParams * getParPtr() {return par_;}

	/** get message pointer */
	virtual StoMessage * getMessagePtr() {return message_;}

protected:

	/** update time stamp and time remains */
	virtual STO_RTN_CODE ticToc();

protected:

	DecModel * model_;     /**< DecModel object */
	DspParams * par_;      /**< parameters */
	StoMessage * message_; /**< message */

	STO_RTN_CODE status_; /**< solution status */
	double * primsol_;    /**< primal solution */
	double * dualsol_;    /**< dual solution */
	double primobj_;      /**< primal objective */
	double dualobj_;      /**< dual objective */

	double time_remains_; /**< time limit */
	double tic_;          /**< time stamp */

public:

	/** solver statistics */
	vector<STO_RTN_CODE> s_statuses_;  /**< history of solution statuses */
	vector<double>       s_primobjs_;  /**< history of primal objective values */
	vector<double>       s_dualobjs_;  /**< history of dual objective values */
	vector<double*>      s_primsols_;  /**< history of primal solutions */
	vector<double>       s_cputimes_;  /**< history of cpu times */
	vector<double>       s_walltimes_; /**< history of wall times */
};

#endif /* DECSOLVER_H_ */
