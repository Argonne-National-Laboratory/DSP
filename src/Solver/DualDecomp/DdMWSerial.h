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
	DdMWSerial();

	/** constructor */
	DdMWSerial(
			DdMaster *        master, /**< master problem */
			vector<DdWorker*> worker  /**< worker for finding lower bounds */);

	/** destructor */
	virtual ~DdMWSerial();

	/** run the framework */
	virtual DSP_RTN_CODE run();

protected:

	/** initialize */
	virtual DSP_RTN_CODE init();

	/** run master process */
	virtual DSP_RTN_CODE runMaster();

	/** run worker processes */
	virtual DSP_RTN_CODE runWorker();

	/** finalize */
	virtual DSP_RTN_CODE finalize();
};

#endif /* SRC_SOLVER_DUALDECOMP_DDMWSERIAL_H_ */
