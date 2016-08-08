/*
 * BaseMasterWorker.h
 *
 *  Created on: Feb 9, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_BASEMASTERWORKER_H_
#define SRC_SOLVER_BASEMASTERWORKER_H_

#include "CoinFinite.hpp"
/** DSP headers */
#include "Utility/DspMacros.h"
#include "Utility/DspRtnCodes.h"
#include "Utility/DspTypes.h"

/**
 * This abstract class defines a basic master-worker framework.
 */
class BaseMasterWorker {
public:

	/** constructor */
	BaseMasterWorker():
		status_(DSP_STAT_UNKNOWN),
		primsol_(NULL),
		dualsol_(NULL),
		primobj_(+COIN_DBL_MAX),
		dualobj_(-COIN_DBL_MAX)
	{}

	/** destructor */
	virtual ~BaseMasterWorker()
	{
		FREE_ARRAY_PTR(primsol_);
		FREE_ARRAY_PTR(dualsol_);
	}

	/** initialize */
	virtual DSP_RTN_CODE init() = 0;

	/** run the framework */
	virtual DSP_RTN_CODE run() = 0;

	/** finalize */
	virtual DSP_RTN_CODE finalize() = 0;

protected:

	/** run master process */
	virtual DSP_RTN_CODE runMaster() = 0;

	/** run worker processes */
	virtual DSP_RTN_CODE runWorker() = 0;

protected:

	/** solution info */
	DSP_RTN_CODE status_; /**< solution status */
	double * primsol_;    /**< primal solution in extensive form */
	double * dualsol_;    /**< dual solution in extensive form */
	double   primobj_;    /**< primal objective bound */
	double   dualobj_;    /**< dual objective bound */
};

#endif /* SRC_SOLVER_BASEMASTERWORKER_H_ */
