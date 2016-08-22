/*
 * DdDriver.h
 *
 *  Created on: Feb 16, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DUALDECOMP_DDDRIVER_H_
#define SRC_SOLVER_DUALDECOMP_DDDRIVER_H_

#include "Solver/DspDriver.h"
#include "Solver/DualDecomp/DdMW.h"

class DdDriver: public DspDriver {
public:

	/** constructor */
	DdDriver(
			DspParams * par,
			DecModel * model);

	/** destructor */
	virtual ~DdDriver();

	/** initialize */
	virtual DSP_RTN_CODE init() {return DSP_RTN_OK;}

	/** run */
	virtual DSP_RTN_CODE run() {return DSP_RTN_OK;}

	/** finalize */
	virtual DSP_RTN_CODE finalize() {return DSP_RTN_OK;}

protected:

	DdMW * mw_;
};

#endif /* SRC_SOLVER_DUALDECOMP_DDDRIVER_H_ */
