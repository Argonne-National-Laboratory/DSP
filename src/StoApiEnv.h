/*
 * StoApiEnv.h
 *
 *  Created on: Nov 6, 2014
 *      Author: kibaekkim
 */

#ifndef STOAPIENV_H_
#define STOAPIENV_H_

#include "Model/TssModel.h"
#include "Solver/StoParam.h"
#include "Solver/TssSolver.h"
#include "Utility/StoMacros.h"

class StoApiEnv
{
public:
	StoApiEnv();
	virtual ~StoApiEnv();

public:
	StoParam * par_;
	TssModel * tss_;
	TssSolver * solver_;
};

#endif /* STOAPIENV_H_ */
