/*
 * DdMaster.cpp
 *
 *  Created on: Feb 9, 2016
 *      Author: kibaekkim
 */

#include <Solver/DualDecomp/DdMaster.h>

DdMaster::DdMaster(
		DspParams *  par,    /**< parameter pointer */
		DecModel *   model,  /**< model pointer */
		DspMessage * message /**< message pointer */):
		//int nworkers          /**< number of workers */):
DdSolver(par, model, message),
si_(NULL),
bestprimobj_(COIN_DBL_MAX),
bestdualobj_(-COIN_DBL_MAX),
bestprimsol_(NULL),
bestdualsol_(NULL),
lambda_(NULL),
subprimobj_(NULL),
subdualobj_(NULL),
subsolution_(NULL)
{}

DdMaster::~DdMaster()
{
	FREE_PTR(si_);
//	FREE_ARRAY_PTR(subindex_);
	FREE_ARRAY_PTR(bestprimsol_);
	FREE_ARRAY_PTR(bestdualsol_);
	FREE_ARRAY_PTR(subprimobj_);
	FREE_ARRAY_PTR(subdualobj_);
	FREE_2D_ARRAY_PTR(model_->getNumSubproblems(),subsolution_);
}

DSP_RTN_CODE DdMaster::init()
{
	BGN_TRY_CATCH

	/** status */
	status_ = DSP_STAT_MW_CONTINUE;

	/** time stamp */
	ticToc();

	/** allocate memory */
//	subindex_    = new int [model_->getNumSubproblems()];
	bestprimsol_ = new double [model_->getFullModelNumCols()];
	bestdualsol_ = new double [model_->getNumCouplingRows()];
	subprimobj_  = new double [model_->getNumSubproblems()];
	subdualobj_  = new double [model_->getNumSubproblems()];
	subsolution_ = new double * [model_->getNumSubproblems()];
	for (int s = 0; s < model_->getNumSubproblems(); ++s)
		subsolution_[s] = new double [model_->getNumSubproblemCouplingCols(s)];

	/** initialize */
	CoinZeroN(bestprimsol_, model_->getFullModelNumCols());
	CoinZeroN(bestdualsol_, model_->getNumCouplingRows());
	CoinZeroN(subprimobj_, model_->getNumSubproblems());
	CoinZeroN(subdualobj_, model_->getNumSubproblems());

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
