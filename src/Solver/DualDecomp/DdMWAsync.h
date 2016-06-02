/*
 * DdMWAsync.h
 *
 *  Created on: Mar 23, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DUALDECOMP_DDMWASYNC_H_
#define SRC_SOLVER_DUALDECOMP_DDMWASYNC_H_

#include <Solver/DualDecomp/DdMWPara.h>

class DdMWAsync: public DdMWPara {

	typedef vector<CoinPackedVector*> Solutions;

public:

	/** constructor */
	DdMWAsync(
			MPI_Comm          comm,   /**< MPI communicator */
			DdMaster *        master, /**< master problem */
			vector<DdWorker*> worker  /**< worker for finding lower bounds */);

	/** destructor */
	virtual ~DdMWAsync();

protected:

	/** initialize */
	virtual DSP_RTN_CODE init();

	/** run master process */
	virtual DSP_RTN_CODE runMaster();

	/** run worker processes */
	virtual DSP_RTN_CODE runWorker();

private:

	/** run master initialization */
	DSP_RTN_CODE runMasterInit();

	/** run master core method */
	DSP_RTN_CODE runMasterCore();

	/** run worker initialization */
	DSP_RTN_CODE runWorkerInit();

	/** run worker core method */
	DSP_RTN_CODE runWorkerCore();

	/** run worker extended methods */
	DSP_RTN_CODE runWorkerExt();

	/** run worker UB method */
	DSP_RTN_CODE runWorkerUB(DdWorkerUB * workerub);

	/** send coupling solution */
	DSP_RTN_CODE sendCouplingSolutions(
			int & maxnum_of_solutions,
			int & maxnum_of_indices,
			int * number_of_elements,
			int * indices,
			double * elements,
			MPI_Request * request);
};

#endif /* SRC_SOLVER_DUALDECOMP_DDMWASYNC_H_ */
