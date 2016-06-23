/*
 * DdMaster.cpp
 *
 *  Created on: Feb 9, 2016
 *      Author: kibaekkim
 */

#include <Solver/DualDecomp/DdMaster.h>

DdMaster::DdMaster(
		DspParams *  par,     /**< parameter pointer */
		DecModel *   model,   /**< model pointer */
		DspMessage * message, /**< message pointer */
		int nworkers          /**< number of workers */):
DdSolver(par, model, message),
si_(NULL),
bestprimobj_(COIN_DBL_MAX),
bestdualobj_(-COIN_DBL_MAX),
nworkers_(nworkers) {}

DdMaster::~DdMaster()
{
	FREE_PTR(si_);
}

DSP_RTN_CODE DdMaster::init()
{
	BGN_TRY_CATCH

	/** status */
	status_ = DSP_STAT_MW_CONTINUE;

	/** time stamp */
	ticToc();

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}


/** set init solution */
DSP_RTN_CODE DdMaster::setInitSolution(const double * sol)
{
	BGN_TRY_CATCH

	if (primsol_)
		CoinCopyN(sol, si_->getNumCols(), primsol_);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}
