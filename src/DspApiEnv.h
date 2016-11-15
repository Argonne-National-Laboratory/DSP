/*
 * StoApiEnv.h
 *
 *  Created on: Nov 6, 2014
 *      Author: kibaekkim
 */

#ifndef DSPAPIENV_H_
#define DSPAPIENV_H_

#include "DspDriver.h"
#include "DecModel.h"
#include "DspParams.h"

class DspApiEnv
{
public:
	DspApiEnv();
	virtual ~DspApiEnv();

public:
	DspDriver * solver_;
	DecModel * model_;
	DspParams * par_;
};

#endif /* DSPAPIENV_H_ */
