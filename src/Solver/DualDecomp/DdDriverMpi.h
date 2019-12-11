/*
 * DdDriverMpi.h
 *
 *  Created on: Jul 8, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DUALDECOMP_DDDRIVERMPI_H_
#define SRC_SOLVER_DUALDECOMP_DDDRIVERMPI_H_

#include "DdDriver.h"

/** A driver class for parallel dual decomposition */
class DdDriverMpi: public DdDriver {
public:

	/** A default constructor. */
	DdDriverMpi(
			DecModel*   model,   /**< model pointer */
			DspParams*  par,     /**< parameters */
			DspMessage* message, /**< message pointer */
			MPI_Comm    comm     /**< MPI communicator */);

	/** A copy constructor. */
	DdDriverMpi(const DdDriverMpi& rhs);

	/** A default destructor. */
	virtual ~DdDriverMpi();

	/** A clone function. */
	virtual DdDriverMpi* clone() const {
		return new DdDriverMpi(*this);
	}

	/** A virtual member for initializing the driver. */
	virtual DSP_RTN_CODE init();

	/** A virtual member for running the driver. */
	virtual DSP_RTN_CODE run();

	/** A virtual member for finalizing the driver. */
	virtual DSP_RTN_CODE finalize();

protected:

	MPI_Comm comm_;
	int comm_rank_;
	int comm_size_;
};

#endif /* SRC_SOLVER_DUALDECOMP_DDDRIVERMPI_H_ */
