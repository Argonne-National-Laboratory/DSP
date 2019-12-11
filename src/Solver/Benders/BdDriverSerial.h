/*
 * BdDriverSerial.h
 *
 *  Created on: Jul 10, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_BENDERS_BDDRIVERSERIAL_H_
#define SRC_SOLVER_BENDERS_BDDRIVERSERIAL_H_

#include "Solver/Benders/BdDriver.h"

/** A driver class for serial Benders decomposition */
class BdDriverSerial: public BdDriver {
public:

	/** A default constructor. */
	BdDriverSerial(
			DecModel *   model,  /**< model pointer */
			DspParams *  par,    /**< parameters */
			DspMessage * message /**< message pointer */);

	/** A copy constructor. */
	BdDriverSerial(const BdDriverSerial& rhs);

	/** A default destructor. */
	virtual ~BdDriverSerial();

	/** A clone function. */
	virtual BdDriverSerial* clone() const {
		return new BdDriverSerial(*this);
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
};

#endif /* SRC_SOLVER_BENDERS_BDDRIVERSERIAL_H_ */
