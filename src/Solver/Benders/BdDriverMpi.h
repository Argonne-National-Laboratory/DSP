/*
 * BdDriverMpi.h
 *
 *  Created on: Jul 10, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_BENDERS_BDDRIVERMPI_H_
#define SRC_SOLVER_BENDERS_BDDRIVERMPI_H_

#include "BdDriver.h"
#include "Solver/Benders/BdMWMpi.h"

/** A driver class for parallel Benders decomposition */
class BdDriverMpi: public BdDriver {
public:

	/** A default constructor. */
	BdDriverMpi(
			DecModel *   model,   /**< model pointer */
			DspParams *  par,     /**< parameters */
			DspMessage * message, /**< message pointer */
			MPI_Comm comm         /**< communicator */);

	/** A copy constructor. */
	BdDriverMpi(const BdDriverMpi& rhs);

	/** A default destructor. */
	virtual ~BdDriverMpi();

	/** A clone function. */
	virtual BdDriverMpi* clone() const {
		return new BdDriverMpi(*this);
	}

	/** A virtual member for initializing the driver. */
	virtual DSP_RTN_CODE init();

	/** A virtual member for running the driver. */
	virtual DSP_RTN_CODE run();

	/** A virtual member for finalizing the driver. */
	virtual DSP_RTN_CODE finalize();

protected:

	/** find lower bound */
	virtual DSP_RTN_CODE findLowerBound();

	/** collect second-stage solutions */
	virtual DSP_RTN_CODE collectSolution();

	MPI_Comm comm_;
	int comm_rank_;
	int comm_size_;
};

#endif /* SRC_SOLVER_BENDERS_BDDRIVERMPI_H_ */
