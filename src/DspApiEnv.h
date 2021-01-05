/*
 * StoApiEnv.h
 *
 *  Created on: Nov 6, 2014
 *      Author: kibaekkim
 */

#ifndef DSPAPIENV_H_
#define DSPAPIENV_H_

#include "DspConfig.h"
#include "Solver/DecSolver.h"
#include "Model/DecModel.h"
#include "Utility/DspParams.h"
#include "Utility/DspMessage.h"

/**
 * A class for DSP API environment.
 */
class DspApiEnv
{
public:
	/** A default constructore */
	DspApiEnv();

	/** A default destructore */
	virtual ~DspApiEnv();

	/** Query DSP version major */
	int getVersionMajor() { return DSP_VERSION_MAJOR; }

	/** Query DSP version minor */
	int getVersionMinor() { return DSP_VERSION_MINOR; }

	/** Query DSP version patch */
	int getVersionPatch() { return DSP_VERSION_PATCH; }

	DecSolver * solver_;   /**< A decomposition solver object */
	DecModel * model_;     /**< A decomposition model object */
	DspParams * par_;      /**< A parameters object */
	DspMessage * message_; /**< A message object */
};

#endif /* DSPAPIENV_H_ */
