/*
 * DecSolver.cpp
 *
 *  Created on: Sep 24, 2014
 *      Author: kibaekkim, ctjandra
 */

#include "DecSolver.h"

DecSolver::DecSolver() :
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

DecSolver::~DecSolver()
{
	FREE_PTR(handler_);
	/** TODO should not release? */
	//FREE_PTR(message_);
	FREE_ARRAY_PTR(solution_);
}

/** load model object; will have a shallow pointer (not deep copy) */
STO_RTN_CODE DecSolver::loadModel(StoParam * par, DecModel * model)
{
	BGN_TRY_CATCH

	if (par == NULL) throw "Error: there is no parameter object.";
	if (model == NULL) throw "Error: there is no model object.";

	par_   = par;
	model_ = model;

	/** number of columns in extensive form */
	int ncols = model_->getFullModelNumCols();

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
