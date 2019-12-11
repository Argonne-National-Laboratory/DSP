/*
 * StoApiEnv.h
 *
 *  Created on: Nov 6, 2014
 *      Author: kibaekkim
 */

#ifndef DSPAPIENV_H_
#define DSPAPIENV_H_

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

	DecSolver * solver_;   /**< A decomposition solver object */
	DecModel * model_;     /**< A decomposition model object */
	DspParams * par_;      /**< A parameters object */
	DspMessage * message_; /**< A message object */
};

#endif /* DSPAPIENV_H_ */
