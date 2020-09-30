/*
 * BdMWMpi.h
 *
 *  Created on: Jul 10, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_BENDERS_BDMWMPI_H_
#define SRC_SOLVER_BENDERS_BDMWMPI_H_

#include "BdMW.h"

class BdMWMpi: public BdMW {
public:

	/** constructor */
	BdMWMpi(
			MPI_Comm comm,
			DecModel *   model,  /**< model pointer */
			DspParams *  par,    /**< parameters */
			DspMessage * message /**< message pointer */);

	/** destructor */
	virtual ~BdMWMpi();

	/** initialize */
	virtual DSP_RTN_CODE init();

	/** run the framework */
	virtual DSP_RTN_CODE run();

	/** finalize */
	virtual DSP_RTN_CODE finalize();

protected:

	/** constraint handler */
	virtual SCIPconshdlrBenders *constraintHandler(bool add_integer_benders);

	/** run master process */
	virtual DSP_RTN_CODE runMaster();

	/** run worker processes */
	virtual DSP_RTN_CODE runWorker();

private:

	/** Common member variables */
	MPI_Comm comm_;
	int comm_rank_;
	int comm_size_;
};

#endif /* SRC_SOLVER_BENDERS_BDMWMPI_H_ */
