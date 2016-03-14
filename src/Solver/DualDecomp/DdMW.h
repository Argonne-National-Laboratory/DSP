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

/**
 * This defines a master-worker framework for dual decomposition.
 */
class DdMW: public BaseMasterWorker {
public:

	/** constructor */
	DdMW(
			MPI_Comm comm,
			DdMaster * master,
			DdWorker * worker);

	/** destructor */
	virtual ~DdMW();

	/** run the framework */
	virtual STO_RTN_CODE run();

protected:

	/** initialize */
	virtual STO_RTN_CODE init();

	/** run master process */
	virtual STO_RTN_CODE runMaster();

	/** run worker processes */
	virtual STO_RTN_CODE runWorker();

	/** run Benders worker processes;
	 * this generates feasiblity/optimality cuts and evaluates upper bound. */
	virtual STO_RTN_CODE runBendersWorker();

	/** finalize */
	virtual STO_RTN_CODE finalize();

private:

	/** run master process (SYNC) */
	virtual STO_RTN_CODE runMasterSync();

	/** run master process (ASYNC) */
	virtual STO_RTN_CODE runMasterAsync();

	/** run worker processes (SYNC) */
	virtual STO_RTN_CODE runWorkerSync();

	/** run worker processes (ASYNC) */
	virtual STO_RTN_CODE runWorkerAsync();

protected:

	DdMaster * master_; /**< master */
	DdWorker * worker_; /**< worker */

private:

	int * nsubprobs_;       /**< number of subproblems for each process */
	int * subprob_indices_; /**< subproblem indices for each process */
	int * subprob_displs_;  /**< displacement of subproblem indices */

	int iteration_limit_;   /**< limit on number of iterations */
};

#endif /* SRC_SOLVER_DECDDMW_H_ */
