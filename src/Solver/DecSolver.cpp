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
	parLogLevel_(0),
	parNodeLim_(COIN_INT_MAX),
	parIterLim_(COIN_INT_MAX),
	parWallLim_(COIN_DBL_MAX),
	parScipDispFreq_(0),
	parScipGapTol_(0),
	parScipTimeLim_(COIN_DBL_MAX),
	parRelaxIntegrality_(NULL),
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
	parRelaxIntegrality_ = NULL;
}

/** load model object; will have a shallow pointer (not deep copy) */
STO_RTN_CODE DecSolver::loadModel(DspParams * par, DecModel * model)
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

	/** parameters */
	parLogLevel_     = par_->getIntParam("LOG_LEVEL");
	parNodeLim_      = par_->getIntParam("NODE_LIM");
	parIterLim_      = par_->getIntParam("ITER_LIM");
	parWallLim_      = par_->getDblParam("WALL_LIM");
	parScipDispFreq_ = par_->getIntParam("SCIP/DISPLAY_FREQ");
	parScipGapTol_   = par_->getDblParam("SCIP/GAP_TOL");
	parScipTimeLim_  = par_->getDblParam("SCIP/TIME_LIM");
	parRelaxIntegrality_ = par_->getBoolPtrParam("RELAX_INTEGRALITY");

	/** create message */
	message_ = new StoMessage(parLogLevel_);

	/** set time limit */
	time_remains_ = parWallLim_;

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
}
