/*
 * DdMWSerial.h
 *
 *  Created on: Apr 8, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DUALDECOMP_DDMWSERIAL_H_
#define SRC_SOLVER_DUALDECOMP_DDMWSERIAL_H_

#include <Solver/DualDecomp/DdMW.h>

class DdMWSerial: public DdMW {
public:

	/** constructor */
	DdMWSerial(
			DecModel *   model,  /**< model pointer */
			DspParams *  par,    /**< parameters */
			DspMessage * message /**< message pointer */);

	/** destructor */
	virtual ~DdMWSerial();

	/** initialize */
	virtual DSP_RTN_CODE init();

	/** run the framework */
	virtual DSP_RTN_CODE run();

	/** finalize */
	virtual DSP_RTN_CODE finalize();

protected:

	/** generate Benders cuts */
	DSP_RTN_CODE generateBendersCuts(
			DdWorkerCGBd * workercg, /**< CG worker pointer */
			Solutions solutions, /**< solutions at which cuts are generated */
			OsiCuts & cuts       /**< cuts generated */);

};

#endif /* SRC_SOLVER_DUALDECOMP_DDMWSERIAL_H_ */
