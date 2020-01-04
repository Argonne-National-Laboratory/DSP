/*
 * BaseMasterWorker.h
 *
 *  Created on: Feb 9, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_BASEMASTERWORKER_H_
#define SRC_SOLVER_BASEMASTERWORKER_H_

#include "CoinFinite.hpp"
#include "Utility/DspMacros.h"
#include "Utility/DspRtnCodes.h"
#include "Utility/DspTypes.h"

/**
 * This abstract class defines a basic master-worker framework.
 */
class BaseMasterWorker {
public:

	/** A default constructor. */
	BaseMasterWorker() {}

	/** A copy constructor. */
	BaseMasterWorker(const BaseMasterWorker& rhs) {}

	/** A default destructor. */
	virtual ~BaseMasterWorker() {}

	/** A clone function */
	virtual BaseMasterWorker* clone() const = 0;

	/** A pure virtual member for initializing the framework. */
	virtual DSP_RTN_CODE init() = 0;

	/** A pure virtual member for running the framework. */
	virtual DSP_RTN_CODE run() = 0;

	/** A pure virtual memeber for finalizing the framework. */
	virtual DSP_RTN_CODE finalize() = 0;

protected:

	/** A pure virtual member to run master process */
	virtual DSP_RTN_CODE runMaster() = 0;

	/** A pure virtual member to run worker processes */
	virtual DSP_RTN_CODE runWorker() = 0;
};

#endif /* SRC_SOLVER_BASEMASTERWORKER_H_ */
