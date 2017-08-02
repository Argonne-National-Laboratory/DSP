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
#include "Solver/DualDecomp/DdWorkerLB.h"
#include "Solver/DualDecomp/DdWorkerUB.h"

#ifndef NO_SCIP
#include "Solver/DualDecomp/DdWorkerCGBd.h"
#endif

/**
 * This defines a master-worker framework for dual decomposition.
 */
class DdMW: public BaseMasterWorker {

protected:

public:

	/** constructor */
	DdMW(
			DecModel *   model,  /**< model pointer */
			DspParams *  par,    /**< parameters */
			DspMessage * message /**< message pointer */);

	/** destructor */
	virtual ~DdMW();

	/** initialize */
	virtual DSP_RTN_CODE init();

	/** run the framework */
	virtual DSP_RTN_CODE run();

	/** finalize */
	virtual DSP_RTN_CODE finalize();

	/** get remaining time */
	virtual double remainingTime() {return parTimeLimit_ - (CoinGetTimeOfDay() - iterstime_);}

protected:

	/** run master process */
	virtual DSP_RTN_CODE runMaster() {return DSP_RTN_OK;}

	/** run worker processes */
	virtual DSP_RTN_CODE runWorker() {return DSP_RTN_OK;}

	/** store coupling solution */
	DSP_RTN_CODE storeCouplingSolutions(Solutions & stored);

	/** check whether solution is duplicate or not; return NULL if duplicate */
	CoinPackedVector * duplicateSolution(
			int size,           /**< size of array */
			const double * x,   /**< current solution */
			Solutions solutions /**< solution pool to check duplication */);

	/** print header info */
	virtual void printHeaderInfo();

	/** print iteration info */
	virtual void printIterInfo();

	/** write output */
	virtual void writeIterInfo(const char * filename);

public:

	DdMaster *        master_; /**< master */
	vector<DdWorker*> worker_; /**< worker for lower bounds */

protected:

	DecModel * model_;     /**< DecModel object */
	DspParams * par_;      /**< parameters */
	DspMessage * message_; /**< message */

	Solutions ubSolutions_; /**< saved solutions that were evaluated for upper bounds */

	OsiCuts * cutsToAdd_; /**< cuts to add */

protected:

	/** parameters */
	int parFeasCuts_;     /**< Benders feasibility cuts */
	int parOptCuts_;      /**< Benders optimality cuts */
	int parEvalUb_;       /**< upper bounds */
	double parTimeLimit_; /**< time limit */

	/** iteration info */
	char   itercode_;
	int    itercnt_;
	double iterstime_; /** start time */

protected:

	vector<double> s_itertime_;    /**< per-iteration time */
	vector<double> s_masterobj_;   /**< master objective */
	vector<double> s_bestprimobj_; /**< best primal objective */
	vector<double> s_bestdualobj_; /**< best dual objective */
};

#endif /* SRC_SOLVER_DECDDMW_H_ */
