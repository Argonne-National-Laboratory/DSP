/*
 * BaseMasterWorker.h
 *
 *  Created on: Feb 9, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_BASEMASTERWORKER_H_
#define SRC_SOLVER_BASEMASTERWORKER_H_

/** DSP headers */
#include <Utility/DspMacros.h>
#include <Utility/DspRtnCodes.h>

/**
 * This abstract class defines a basic master-worker framework.
 */
class BaseMasterWorker {
public:

	/** constructor */
	BaseMasterWorker() {}

	/** destructor */
	virtual ~BaseMasterWorker() {}

	/** run the framework */
	virtual DSP_RTN_CODE run() = 0;

protected:

	/** initialize */
	virtual DSP_RTN_CODE init() = 0;

	/** run master process */
	virtual DSP_RTN_CODE runMaster() = 0;

	/** run worker processes */
	virtual DSP_RTN_CODE runWorker() = 0;

	/** finalize */
	virtual DSP_RTN_CODE finalize() = 0;
};

#endif /* SRC_SOLVER_BASEMASTERWORKER_H_ */