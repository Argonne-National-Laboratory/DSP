/*
 * DecSolver.cpp
 *
 *  Created on: Sep 24, 2014
 *      Author: kibaekkim, ctjandra
 */

#include "Solver/DecSolver.h"

DecSolver::DecSolver(DspParams * par, DecModel * model, StoMessage * message):
	model_(model), par_(par), message_(message),
	status_(STO_STAT_UNKNOWN), primsol_(NULL), dualsol_(NULL),
	primobj_(COIN_DBL_MAX), dualobj_(-COIN_DBL_MAX),
	time_remains_(COIN_DBL_MAX), tic_(0.0) {}

DecSolver::~DecSolver()
{
	FREE_ARRAY_PTR(primsol_);
	for (unsigned i = 0; i < s_primsols_.size(); ++i)
		FREE_ARRAY_PTR(s_primsols_[i]);
	message_ = NULL;
	par_ = NULL;
	model_ = NULL;
}

/** update time stamp and time remains */
STO_RTN_CODE DecSolver::ticToc()
{
	BGN_TRY_CATCH

	time_remains_ -= CoinGetTimeOfDay() - tic_;
	tic_ = CoinGetTimeOfDay();

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
}
