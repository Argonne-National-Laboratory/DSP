/*
 * DdDriverSerial.h
 *
 *  Created on: Jul 8, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DUALDECOMP_DDDRIVERSERIAL_H_
#define SRC_SOLVER_DUALDECOMP_DDDRIVERSERIAL_H_

#include "DdDriver.h"

/** A driver class for serial dual decomposition */
class DdDriverSerial: public DdDriver {
public:

	/** A default constructor. */
	DdDriverSerial(
			DecModel *   model,  /**< model pointer */
			DspParams *  par,    /**< parameters */
			DspMessage * message /**< message pointer */);
	
	/** A copy constructor. */
	DdDriverSerial(const DdDriverSerial& rhs);

	/** A default destructor. */
	virtual ~DdDriverSerial();

	/** A clone function */
	virtual DdDriverSerial* clone() const {
		return new DdDriverSerial(*this);
	}

	/** A virtual member for initializing the driver. */
	virtual DSP_RTN_CODE init();

	/** A virtual member for running the driver. */
	virtual DSP_RTN_CODE run();

	/** A virtual memeber for finalizing the driver. */
	virtual DSP_RTN_CODE finalize();

};

#endif /* SRC_SOLVER_DUALDECOMP_DDDRIVERSERIAL_H_ */
