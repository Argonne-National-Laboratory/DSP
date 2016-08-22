/*
 * BdDriverSerial.h
 *
 *  Created on: Jul 10, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_BENDERS_BDDRIVERSERIAL_H_
#define SRC_SOLVER_BENDERS_BDDRIVERSERIAL_H_

#include "Solver/Benders/BdDriver.h"

class BdDriverSerial: public BdDriver {
public:

	/** constructor */
	BdDriverSerial(
			DspParams * par, /**< parameter pointer */
			DecModel * model /**< model pointer */);

	/** destructor */
	virtual ~BdDriverSerial();

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
};

#endif /* SRC_SOLVER_BENDERS_BDDRIVERSERIAL_H_ */
