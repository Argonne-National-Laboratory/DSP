/*
 * StoApiEnv.h
 *
 *  Created on: Nov 6, 2014
 *      Author: kibaekkim
 */

#ifndef DSPAPIENV_H_
#define DSPAPIENV_H_

#include "Solver/DecSolver.h"
#include "DecModel.h"
#include "DspParams.h"
#include "Utility/DspMessage.h"

class DspApiEnv
{
public:
	DspApiEnv();
	virtual ~DspApiEnv();

public:
	DecSolver * solver_;
	DecModel * model_;
	DspParams * par_;
	DspMessage * message_;
};

#endif /* DSPAPIENV_H_ */
