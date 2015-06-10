/*
 * TssSolver.cpp
 *
 *  Created on: Sep 24, 2014
 *      Author: kibaekkim
 */

#include "TssSolver.h"

TssSolver::TssSolver() :
	model_(NULL),
	par_(NULL),
	message_(NULL),
	time_remains_(COIN_DBL_MAX),
	status_(STO_STAT_UNKNOWN),
	solution_(NULL),
	primalBound_(COIN_DBL_MAX),
	dualBound_(-COIN_DBL_MAX),
	numIterations_(0),
	numNodes_(0),
	clockType_(1),
	solutionTime_(0.)
{
	handler_ = new CoinMessageHandler;
}

TssSolver::~TssSolver()
{
	FREE_PTR(handler_);
	/** TODO should not release? */
	//FREE_PTR(message_);
	FREE_ARRAY_PTR(solution_);
}

/** load model object; will have a shallow pointer (not deep copy) */
STO_RTN_CODE TssSolver::loadModel(StoParam * par, TssModel * model)
{
	BGN_TRY_CATCH

	if (par == NULL) throw "Error: there is no parameter object.";
	if (model == NULL) throw "Error: there is no model object.";

	par_   = par;
	model_ = model;

	/** number of columns in extensive form */
	int ncols = model_->getNumCols(0) + model_->getNumCols(1) * model_->getNumScenarios();

	/** allocate memory for solution */
	assert(solution_ == NULL);
	solution_ = new double [ncols];

	/** create message */
	message_ = new StoMessage;

	/** set time limit */
	time_remains_ = par_->wtimeLimit_;

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
}

/**
 * This converts ctype to the number of integer variables and their indices.
 */
void TssSolver::convertColTypes(int ncols, char * ctype, int & len, int *& ind)
{
	len = 0;

	/** count integer variables */
	for (int j = 0; j < ncols; ++j)
	{
#if RELAX_FIRST_STAGE_INTEGER == 1
		if (j < model_->getNumCols(0)) continue;
#endif
#if RELAX_SECOND_STAGE_INTEGER == 1
		if (j >= model_->getNumCols(0) && j < ncols) continue;
#endif
		if (ctype[j] != 'C')
		{
			len++;
		}
	}

	/** get integer variable indices */
	ind = new int [len];
	for (int j = 0, jj = 0; j < ncols; ++j)
	{
#if RELAX_FIRST_STAGE_INTEGER == 1
		if (j < model_->getNumCols(0)) continue;
#endif
#if RELAX_SECOND_STAGE_INTEGER == 1
		if (j >= model_->getNumCols(0) && j < ncols) continue;
#endif
		if (ctype[j] != 'C')
		{
			ind[jj++] = j;
		}
	}
}

