/*
 * DdDriver.h
 *
 *  Created on: Feb 16, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DUALDECOMP_DDDRIVER_H_
#define SRC_SOLVER_DUALDECOMP_DDDRIVER_H_

#include "Solver/DspDriver.h"
#include "Solver/DualDecomp/DdMWSync.h"
#include "Solver/DualDecomp/DdMWAsync.h"

class DdDriver: public DspDriver {
public:

	/** constructor */
	DdDriver(DspParams * par, DecModel * model);

	/** constructor with MPI */
	DdDriver(DspParams * par, DecModel * model, MPI_Comm comm);

	/** destructor */
	virtual ~DdDriver();

	/** initialize */
	virtual DSP_RTN_CODE init();

	/** run */
	virtual DSP_RTN_CODE run();

	/** finalize */
	virtual DSP_RTN_CODE finalize();

public:

	/** get pointer to master */
	const DdMaster * getMasterPtr() {return master_;}

	/** get pointer to worker */
	vector<DdWorker*> getWorkerList() {return worker_;}

	/** get number of infeasible solutions evaluated */
	int getNumInfeasibleSolutions() {return num_infeasible_solutions_;}

private:

	MPI_Comm comm_;
	int comm_rank_;
	int comm_size_;

	DdMW * mw_;
	DdMaster *        master_; /**< master */
	vector<DdWorker*> worker_; /**< worker place holder */
//	DdWorker *   worker_;   /**< worker for lower bounds */
//	DdWorkerUB * workerUB_; /**< worker for upper bounds */

	int num_infeasible_solutions_; /**< number of infeasible solutions found in subproblems */

};

#endif /* SRC_SOLVER_DUALDECOMP_DDDRIVER_H_ */
