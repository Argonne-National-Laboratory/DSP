/*
 * DeDriver.h
 *
 *  Created on: Feb 17, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DETERMINISTIC_DEDRIVER_H_
#define SRC_SOLVER_DETERMINISTIC_DEDRIVER_H_

#include "Solver/DecSolver.h"

/**
 * This class defines a driver for solving a deterministic equivalent problem.
 */
class DeDriver: public DecSolver {
public:

	/** A default constructor. */
	DeDriver(
			DecModel *   model,  /**< model pointer */
			DspParams *  par,    /**< parameters */
			DspMessage * message /**< message pointer */);

	/** A copy constructor. */
	DeDriver(const DeDriver& rhs);

	/** A default destructor. */
	virtual ~DeDriver();

	/** A clone function. */
	virtual DeDriver* clone() const {
		return new DeDriver(*this);
	}

	/** A virtual member for initilizing the driver. */
	virtual DSP_RTN_CODE init();

	/** A virtual member for running the driver. */
	virtual DSP_RTN_CODE run();
	virtual DSP_RTN_CODE solve();

	/** A virtual memeber for finalizing the driver. */
	virtual DSP_RTN_CODE finalize();

	/** write extensive form in MPS */
	virtual void writeExtMps(const char * name);

private:

	virtual DspOsi * createDspOsi();
};

#endif /* SRC_SOLVER_DETERMINISTIC_DEDRIVER_H_ */
