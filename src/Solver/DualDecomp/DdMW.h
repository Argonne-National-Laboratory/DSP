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

	typedef DdMaster Master;
	typedef DdWorker Worker;

public:

	/** constructor */
	DdMW(
			MPI_Comm comm,
			Master * master,
			Worker * worker);

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

protected:

	Master * master_; /**< master */
	Worker * worker_; /**< worker */
};

#endif /* SRC_SOLVER_DECDDMW_H_ */
