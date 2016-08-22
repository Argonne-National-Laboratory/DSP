/*
 * DdDriverSerial.h
 *
 *  Created on: Jul 8, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DUALDECOMP_DDDRIVERSERIAL_H_
#define SRC_SOLVER_DUALDECOMP_DDDRIVERSERIAL_H_

#include "DdDriver.h"

class DdDriverSerial: public DdDriver {
public:

	/** constructor */
	DdDriverSerial(
			DspParams * par,
			DecModel * model);

	/** destructor */
	virtual ~DdDriverSerial() {}

	/** initialize */
	virtual DSP_RTN_CODE init();

	/** run */
	virtual DSP_RTN_CODE run();

	/** finalize */
	virtual DSP_RTN_CODE finalize();

};

#endif /* SRC_SOLVER_DUALDECOMP_DDDRIVERSERIAL_H_ */
