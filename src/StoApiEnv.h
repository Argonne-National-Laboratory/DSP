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
#include "Solver/StoParam.h"
#include "Solver/DecSolver.h"
#include "Utility/StoMacros.h"

class StoApiEnv
{
public:
	StoApiEnv();
	virtual ~StoApiEnv();

public:
	StoParam * par_;
	DecModel * model_;
	DecSolver * solver_;
	DecData * decdata_;  /**< decomposition data: only used when a general decomposable model is supplied */

};

#endif /* STOAPIENV_H_ */
