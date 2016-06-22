/*
 * DdMasterSync.cpp
 *
 *  Created on: Mar 25, 2016
 *      Author: kibaekkim
 */

#include <Solver/DualDecomp/DdMasterSync.h>

DdMasterSync::DdMasterSync(
		DspParams *  par,     /**< parameter pointer */
		DecModel *   model,   /**< model pointer */
		DspMessage * message, /**< message pointer */
		int nworkers          /**< number of workers */):
DdMaster(par, model, message, nworkers),
worker_(-1),
nsubprobs_(0),
subindex_(NULL),
subprimobj_(NULL),
subdualobj_(NULL),
subsolution_(NULL)
{
	// TODO Auto-generated constructor stub

}

DdMasterSync::~DdMasterSync()
{
	FREE_ARRAY_PTR(subindex_);
	FREE_ARRAY_PTR(subprimobj_);
	FREE_ARRAY_PTR(subdualobj_);
	FREE_2D_ARRAY_PTR(model_->getNumSubproblems(),subsolution_);
}

DSP_RTN_CODE DdMasterSync::init()
{
	BGN_TRY_CATCH

	DdMaster::init();

	/** allocate memory */
	subindex_    = new int [model_->getNumSubproblems()];
	subprimobj_  = new double [model_->getNumSubproblems()];
	subdualobj_  = new double [model_->getNumSubproblems()];
	subsolution_ = new double * [model_->getNumSubproblems()];
	for (int s = 0; s < model_->getNumSubproblems(); ++s)
		subsolution_[s] = new double [model_->getNumSubproblemCouplingCols(s)];

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}
