/*
 * BaseMasterWorker.cpp
 *
 *  Created on: Feb 18, 2016
 *      Author: kibaekkim
 */

#include "Solver/BaseMasterWorker.h"


BaseMasterWorker::BaseMasterWorker(MPI_Comm comm):
	comm_(comm), sync_(true),
	scount_(0), rcount_(0), sendbuf_(NULL), recvbuf_(NULL),
	scounts_(NULL), sdispls_(NULL), rcounts_(NULL), rdispls_(NULL)
{
	MPI_Comm_rank(comm, &comm_rank_);
	MPI_Comm_size(comm, &comm_size_);
}

BaseMasterWorker::~BaseMasterWorker()
{
	FREE_ARRAY_PTR(sendbuf_);
	FREE_ARRAY_PTR(recvbuf_);
	FREE_ARRAY_PTR(scounts_);
	FREE_ARRAY_PTR(sdispls_);
	FREE_ARRAY_PTR(rcounts_);
	FREE_ARRAY_PTR(rdispls_);
}

STO_RTN_CODE BaseMasterWorker::run()
{
	BGN_TRY_CATCH

	/** initialize */
	init();

	/** run master process */
	runMaster();

	/** run worker processes */
	runWorker();

	/** finalize */
	finalize();

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
}


