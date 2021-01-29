/*
 * BdMWSerial.h
 *
 *  Created on: Jul 10, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_BENDERS_BDMWSERIAL_H_
#define SRC_SOLVER_BENDERS_BDMWSERIAL_H_

#include "BdMW.h"

class BdMWSerial: public BdMW {
public:

	/** constructor */
	BdMWSerial(
			DecModel *   model,  /**< model pointer */
			DspParams *  par,    /**< parameters */
			DspMessage * message /**< message pointer */);

	/** destructor */
	virtual ~BdMWSerial();

	/** initialize */
	virtual DSP_RTN_CODE init();

	/** run the framework */
	virtual DSP_RTN_CODE run();

	/** finalize */
	virtual DSP_RTN_CODE finalize();

protected:

	/** constraint handler */
	virtual SCIPconshdlrBenders * constraintHandler();

	/** Benders Callback Function */
	virtual BendersCallback * BendersCallbackFunc();
};

#endif /* SRC_SOLVER_BENDERS_BDMWSERIAL_H_ */
