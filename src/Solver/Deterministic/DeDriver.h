/*
 * DeDriver.h
 *
 *  Created on: Feb 17, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DETERMINISTIC_DEDRIVER_H_
#define SRC_SOLVER_DETERMINISTIC_DEDRIVER_H_

#include "Solver/DspDriver.h"
#include "SolverInterface/SolverInterface.h"

/**
 * This class defines a driver for solving a deterministic equivalent problem.
 */
class DeDriver: public DspDriver {
public:

	/** constructor */
	DeDriver(DspParams * par, DecModel * model);

	/** destructor */
	virtual ~DeDriver();

	/** initilize */
	virtual DSP_RTN_CODE init();

	/** run */
	virtual DSP_RTN_CODE run();

	/** finalize */
	virtual DSP_RTN_CODE finalize();

private:

	SolverInterface * si_; /**< my solver interface */
};

#endif /* SRC_SOLVER_DETERMINISTIC_DEDRIVER_H_ */
