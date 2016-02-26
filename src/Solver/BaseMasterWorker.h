/*
 * BaseMasterWorker.h
 *
 *  Created on: Feb 9, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_BASEMASTERWORKER_H_
#define SRC_SOLVER_BASEMASTERWORKER_H_

/** MPI */
#include "mpi.h"

/** DSP headers */
#include "StoRtnCodes.h"
#include "StoMacros.h"

/**
 * This abstract class defines a basic master-worker framework.
 */
class BaseMasterWorker {
public:

	/** constructor */
	BaseMasterWorker(MPI_Comm comm);

	/** destructor */
	virtual ~BaseMasterWorker();

	/** run the framework */
	virtual STO_RTN_CODE run();

protected:

	/** initialize */
	virtual STO_RTN_CODE init() = 0;

	/** run master process */
	virtual STO_RTN_CODE runMaster() = 0;

	/** run worker processes */
	virtual STO_RTN_CODE runWorker() = 0;

	/** finalize */
	virtual STO_RTN_CODE finalize() = 0;

protected:

	/** Common member variables */
	MPI_Comm comm_;
	int comm_rank_;
	int comm_size_;
	bool sync_; /**< indicate whether parallelism is synchronous or asynchronous */

	int scount_; /**< sending message size */
	int rcount_; /**< receiving message size */
	double * sendbuf_; /**< sending message buffer */
	double * recvbuf_; /**< receiving message buffer */

	int * scounts_; /**< sending message size array for MPI_Scatterv */
	int * sdispls_; /**< receiving message displacement for MPI_Scatterv */

	int * rcounts_; /**< receiving message size array for MPI_Gatherv */
	int * rdispls_; /**< receiving message displacement for MPI_Gatherv */

};

#endif /* SRC_SOLVER_BASEMASTERWORKER_H_ */
