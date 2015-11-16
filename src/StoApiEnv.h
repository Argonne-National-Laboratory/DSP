/*
 * StoApiEnv.h
 *
 *  Created on: Nov 6, 2014
 *      Author: kibaekkim
 */

#ifndef STOAPIENV_H_
#define STOAPIENV_H_

#include "Model/TssModel.h"
#include "Model/DecData.h"
#include "Utility/DspParams.h"
#include "Solver/DecSolver.h"
#include "Utility/StoMacros.h"

class StoApiEnv
{
public:
	StoApiEnv();
	virtual ~StoApiEnv();

public:
	DecSolver * solver_;
	DecData * decdata_;  /**< decomposition data: only used when a general decomposable model is supplied */
	DecModel * model_;
	DspParams * par_;

};

#endif /* STOAPIENV_H_ */
