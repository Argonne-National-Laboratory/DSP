/*
 * DdDriverMpi.h
 *
 *  Created on: Jul 8, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DUALDECOMP_DDDRIVERMPI_H_
#define SRC_SOLVER_DUALDECOMP_DDDRIVERMPI_H_

#include "DdDriver.h"

class DdDriverMpi: public DdDriver {
public:

	/** constructor */
	DdDriverMpi(
			DspParams * par,
			DecModel * model,
			MPI_Comm comm);

	/** destructor */
	virtual ~DdDriverMpi() {}

	/** initialize */
	virtual DSP_RTN_CODE init();

	/** run */
	virtual DSP_RTN_CODE run();

	/** finalize */
	virtual DSP_RTN_CODE finalize();

private:

	MPI_Comm comm_;
	int comm_rank_;
	int comm_size_;
};

#endif /* SRC_SOLVER_DUALDECOMP_DDDRIVERMPI_H_ */
