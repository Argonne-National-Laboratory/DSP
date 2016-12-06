/*
 * DspDriverMpi.h
 *
 *  Created on: Dec 5, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_DSPDRIVERMPI_H_
#define SRC_DSPDRIVERMPI_H_

#include "mpi.h"
#include <DspDriver.h>

class DspDriverMpi: public DspDriver {
public:
	DspDriverMpi(
			DecModel* model, /**< model pointer */
			DspParams* par,  /**< parameter pointer */
			MPI_Comm comm    /**< MPI Communicator */);

	virtual ~DspDriverMpi() {
		comm_ = MPI_COMM_NULL;
	}

	/** initialize */
	virtual DSP_RTN_CODE init();

	/** run */
	virtual DSP_RTN_CODE run();

protected:

	MPI_Comm comm_;
	int comm_size_;
	int comm_rank_;
};

#endif /* SRC_DSPDRIVERMPI_H_ */
