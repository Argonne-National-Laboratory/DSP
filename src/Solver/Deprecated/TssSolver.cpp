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
	parLogLevel_(0),
	parNodeLim_(COIN_INT_MAX),
	parIterLim_(COIN_INT_MAX),
	parWallLim_(COIN_DBL_MAX),
	parScipDispFreq_(0),
	parScipGapTol_(0),
	parScipTimeLim_(COIN_DBL_MAX),
	parRelaxIntegrality_(NULL),
	status_(DSP_STAT_UNKNOWN),
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
	parRelaxIntegrality_ = NULL;
}

/** load model object; will have a shallow pointer (not deep copy) */
DSP_RTN_CODE TssSolver::loadModel(DspParams * par, TssModel * model)
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

	/** parameters */
	parLogLevel_         = par_->getIntParam("LOG_LEVEL");
	parNodeLim_          = par_->getIntParam("NODE_LIM");
	parIterLim_          = par_->getIntParam("ITER_LIM");
	parWallLim_          = par_->getDblParam("WALL_LIM");
	parScipDispFreq_     = par_->getIntParam("SCIP/DISPLAY_FREQ");
	parScipGapTol_       = par_->getDblParam("SCIP/GAP_TOL");
	parScipTimeLim_      = par_->getDblParam("SCIP/TIME_LIM");
	parRelaxIntegrality_ = par_->getBoolPtrParam("RELAX_INTEGRALITY");

	/** create message */
	message_ = new DspMessage(parLogLevel_);

	/** set time limit */
	time_remains_ = parWallLim_;

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
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

