/*
 * DwSolverMpi.h
 *
 *  Created on: Dec 5, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DANTZIGWOLFE_DWSOLVERMPI_H_
#define SRC_SOLVER_DANTZIGWOLFE_DWSOLVERMPI_H_

#include "Utility/DspMpi.h"
#include "Solver/DantzigWolfe/DwSolverSerial.h"

class DwSolverMpi: public DwSolverSerial {
public:

    /** default constructor */
	DwSolverMpi(
			DecModel*   model,   /**< model pointer */
			DspParams*  par,     /**< parameters */
			DspMessage* message, /**< message pointer */
			MPI_Comm    comm     /**< MPI communicator */);

	virtual ~DwSolverMpi();

	/** initialize */
	virtual DSP_RTN_CODE init();

	/** solve */
	virtual DSP_RTN_CODE solve();

	/** finalize */
	virtual DSP_RTN_CODE finalize();

protected:

	MPI_Comm comm_;
	int comm_size_;
	int comm_rank_;
};

#endif /* SRC_SOLVER_DANTZIGWOLFE_DWSOLVERMPI_H_ */
