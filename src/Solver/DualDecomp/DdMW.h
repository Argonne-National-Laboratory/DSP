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
#include "Solver/DualDecomp/DdWorkerUBQcp.h"
#include "Solver/DualDecomp/DdDroWorkerUB.h"

#ifdef DSP_HAS_SCIP
#include "Solver/DualDecomp/DdWorkerCGBd.h"
#endif

/** A base master-worker class for dual decomposition */
class DdMW: public BaseMasterWorker {

protected:

public:

	/** A default constructor. */
	DdMW(
			DecModel *   model,  /**< model pointer */
			DspParams *  par,    /**< parameters */
			DspMessage * message /**< message pointer */);

	/** A copy constructor. */
	DdMW(const DdMW& rhs);

	/** A default destructor. */
	virtual ~DdMW();

	/** A clone function */
	virtual DdMW* clone() const {
		return new DdMW(*this);
	}

	/** A virtual member for initializing the framework. */
	virtual DSP_RTN_CODE init();

	/** A virtual member for running the framework. */
	virtual DSP_RTN_CODE run();

	/** A virtual memeber for finalizing the framework. */
	virtual DSP_RTN_CODE finalize();

	/** get number of iterations */
	virtual int getIterationCount() {return itercnt_;}

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

	/** parameters */
	int parFeasCuts_;     /**< Benders feasibility cuts */
	int parOptCuts_;      /**< Benders optimality cuts */
	int parEvalUb_;       /**< upper bounds */
	double parTimeLimit_; /**< time limit */

	/** iteration info */
	char   itercode_;
	int    itercnt_;
	double iterstime_; /** start time */

	vector<double> s_itertime_;    /**< per-iteration time */
	vector<double> s_masterobj_;   /**< master objective */
	vector<double> s_bestprimobj_; /**< best primal objective */
	vector<double> s_bestdualobj_; /**< best dual objective */
};

#endif /* SRC_SOLVER_DECDDMW_H_ */
