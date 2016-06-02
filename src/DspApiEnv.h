/*
 * StoApiEnv.h
 *
 *  Created on: Nov 6, 2014
 *      Author: kibaekkim
 */

#ifndef DSPAPIENV_H_
#define DSPAPIENV_H_

#include <Utility/DspMacros.h>
#include "Model/TssModel.h"
#include "Model/DecData.h"
#include "Utility/DspParams.h"
#include "Solver/DspDriver.h"
#include "Solver/DecSolver.h"

class DspApiEnv
{
public:
	DspApiEnv();
	virtual ~DspApiEnv();

public:
	DspDriver * solver_;
	DecData * decdata_;  /**< decomposition data: only used when a general decomposable model is supplied */
	DecModel * model_;
	DspParams * par_;

};

#endif /* DSPAPIENV_H_ */
