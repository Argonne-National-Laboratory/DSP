/*
 * BdDriverMpi.h
 *
 *  Created on: Jul 10, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_BENDERS_BDDRIVERMPI_H_
#define SRC_SOLVER_BENDERS_BDDRIVERMPI_H_

#include "BdDriver.h"

class BdDriverMpi: public BdDriver {
public:

	/** constructor */
	BdDriverMpi(DspParams * par, DecModel * model, MPI_Comm comm);

	/** destructor */
	virtual ~BdDriverMpi();

	/** initialize */
	virtual DSP_RTN_CODE init();

	/** run */
	virtual DSP_RTN_CODE run();

	/** finalize */
	virtual DSP_RTN_CODE finalize();

protected:

	/** find lower bound */
	virtual DSP_RTN_CODE findLowerBound();

	/** collect second-stage solutions */
	virtual DSP_RTN_CODE collectSolution();

private:

	MPI_Comm comm_;
	int comm_rank_;
	int comm_size_;
};

#endif /* SRC_SOLVER_BENDERS_BDDRIVERMPI_H_ */
