/*
 * DdMasterSync.cpp
 *
 *  Created on: Mar 25, 2016
 *      Author: kibaekkim
 */

#include <Solver/DualDecomp/DdMasterSync.h>

DdMasterSync::DdMasterSync(
		DspParams * par,
		DecModel * model,
		DspMessage * message,
		int nworkers,
		int maxnumsubprobs):
	DdMaster(par, model, message, nworkers, maxnumsubprobs),
	worker_(-1), nsubprobs_(0), subindex_(NULL), subprimobj_(NULL), subdualobj_(NULL), subsolution_(NULL)
{
	// TODO Auto-generated constructor stub

}

DdMasterSync::~DdMasterSync()
{
	FREE_ARRAY_PTR(subindex_);
	FREE_ARRAY_PTR(subprimobj_);
	FREE_ARRAY_PTR(subdualobj_);
	FREE_2D_ARRAY_PTR(maxnumsubprobs_,subsolution_);
}

DSP_RTN_CODE DdMasterSync::init()
{
	BGN_TRY_CATCH

	DdMaster::init();

	/** allocate memory */
	subindex_    = new int [maxnumsubprobs_];
	subprimobj_  = new double [maxnumsubprobs_];
	subdualobj_  = new double [maxnumsubprobs_];
	subsolution_ = new double * [maxnumsubprobs_];
	for (int s = 0; s < maxnumsubprobs_; ++s)
		subsolution_[s] = new double [model_->getNumSubproblemCouplingCols(s)];

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}
