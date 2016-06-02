/*
 * DdMW.h
 *
 *  Created on: Feb 9, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DECDDMW_H_
#define SRC_SOLVER_DECDDMW_H_

#include "Solver/BaseMasterWorker.h"
#include "Solver/DualDecomp/DdMaster.h"
#include "Solver/DualDecomp/DdWorker.h"
#include "Solver/DualDecomp/DdWorkerCGBd.h"
#include "Solver/DualDecomp/DdWorkerLB.h"
#include "Solver/DualDecomp/DdWorkerUB.h"

/**
 * This defines a master-worker framework for dual decomposition.
 */
class DdMW: public BaseMasterWorker {

	typedef vector<CoinPackedVector*> Solutions;

public:

	/** constructor */
	DdMW(
			DdMaster *        master, /**< master problem */
			vector<DdWorker*> worker  /**< worker for finding lower bounds */);

	/** destructor */
	virtual ~DdMW();

	/** initialize */
	virtual DSP_RTN_CODE init();

	/** run the framework */
	virtual DSP_RTN_CODE run();

	/** finalize */
	virtual DSP_RTN_CODE finalize();

protected:

	/** run master process */
	virtual DSP_RTN_CODE runMaster() {return DSP_RTN_OK;}

	/** run worker processes */
	virtual DSP_RTN_CODE runWorker() {return DSP_RTN_OK;}

	/** check whether solution is duplicate or not; return NULL if duplicate */
	CoinPackedVector * duplicateSolution(
			int size,           /**< size of array */
			const double * x,   /**< current solution */
			Solutions solutions /**< solution pool to check duplication */);

	/** print iteration info */
	virtual void printIterInfo();

protected:

	DdMaster *        master_; /**< master */
	vector<DdWorker*> worker_; /**< worker for lower bounds */

	DecModel * model_;     /**< DecModel object */
	DspParams * par_;      /**< parameters */
	DspMessage * message_; /**< message */

	Solutions ubSolutions_; /**< saved solutions that were evaluated for upper bounds */

	OsiCuts * cutsToAdd_; /**< cuts to add */

protected:

	/** parameters */
	int parFeasCuts_; /**< Benders feasibility cuts */
	int parOptCuts_;  /**< Benders optimality cuts */
	int parEvalUb_;   /**< upper bounds */

	/** iteration info */
	char itercode_;
	int itercnt_;
	double iterstime_; /** start time */
};

#endif /* SRC_SOLVER_DECDDMW_H_ */
