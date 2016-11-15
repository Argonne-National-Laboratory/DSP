/*
 * BdMaster.h
 *
 *  Created on: Feb 18, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_BENDERS_BDMASTER_H_
#define SRC_SOLVER_BENDERS_BDMASTER_H_

/** DSP */
#include "Solver/DecSolver.h"
#include "Solver/Benders/BdWorker.h"

class BdMaster: public DecSolver {

public:

	/** constructor */
	BdMaster(DecModel * model, DspParams * par, DspMessage * message);

	/** destructor */
	virtual ~BdMaster();

	/** initialize */
	virtual DSP_RTN_CODE init();

	/** solve */
	virtual DSP_RTN_CODE solve();

	/** finalize */
	virtual DSP_RTN_CODE finalize();

	/** set worker pointer */
	void setWorkerPtr(DdWorker* worker) {worker_ = worker;}

public:

	/** set objective bounds */
	virtual DSP_RTN_CODE setObjectiveBounds(double upper, double lower);

	/** set auxiliary column data */
	virtual DSP_RTN_CODE setAuxVarData(int size, double * obj, double * clbd, double * cubd);

	/** set initial solutions */
	virtual DSP_RTN_CODE setSolutions(Solutions initsols);

protected:

	/** create problem */
	virtual DSP_RTN_CODE createProblem();

protected:

	BdWorker* worker_; /**< Benders worker */

	int      naux_;     /**< number of auxiliary variables */
	double * obj_aux_;  /**< auxiliary variable objectives */
	double * clbd_aux_; /**< auxiliary variable lower bounds */
	double * cubd_aux_; /**< auxiliary variable upper bounds */
};

#endif /* SRC_SOLVER_BENDERS_BDMASTER_H_ */
