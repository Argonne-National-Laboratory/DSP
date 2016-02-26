/*
 * DdMaster.cpp
 *
 *  Created on: Feb 9, 2016
 *      Author: kibaekkim
 */

#include <Solver/DualDecomp/DdMaster.h>

DdMaster::DdMaster(DspParams * par, DecModel * model, StoMessage * message, int nworkers, int maxnumsubprobs) :
	DdSolver(par, model, message),
	si_(NULL), bestprimobj_(COIN_DBL_MAX), bestdualobj_(-COIN_DBL_MAX),
	nworkers_(nworkers), worker_(-1), nsubprobs_(0), maxnumsubprobs_(maxnumsubprobs),
	subindex_(NULL), subprimobj_(NULL), subdualobj_(NULL), subsolution_(NULL) {}

DdMaster::~DdMaster()
{
	FREE_PTR(si_);
	FREE_ARRAY_PTR(subindex_);
	FREE_ARRAY_PTR(subprimobj_);
	FREE_ARRAY_PTR(subdualobj_);
	FREE_2D_ARRAY_PTR(maxnumsubprobs_,subsolution_);
}

STO_RTN_CODE DdMaster::init()
{
	BGN_TRY_CATCH

	/** message count */
	scount_ = 1;
	rcount_ = 1;
	for (unsigned s = 0; s < maxnumsubprobs_; ++s)
	{
		scount_ += 1 + model_->getNumSubproblemCouplingRows(s);
		rcount_ += 3 + model_->getNumSubproblemCouplingCols(s);
	}
	/** message buffer */
	sbuf_ = new double [scount_];
	sbuf_[0] = static_cast<double>(STO_STAT_MW_CONTINUE);

	/** allocate memory */
	subindex_    = new int [maxnumsubprobs_];
	subprimobj_  = new double [maxnumsubprobs_];
	subdualobj_  = new double [maxnumsubprobs_];
	subsolution_ = new double * [maxnumsubprobs_];
	for (int s = 0; s < maxnumsubprobs_; ++s)
		subsolution_[s] = new double [model_->getNumSubproblemCouplingCols(s)];

	/** status */
	status_ = STO_STAT_MW_CONTINUE;

	/** time stamp */
	ticToc();

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
}

STO_RTN_CODE DdMaster::recvMessage(int source, int size, double* message)
{
	BGN_TRY_CATCH

	int pos = 0;
	nsubprobs_ = static_cast<int>(message[pos++]);
	for (int s = 0; s < nsubprobs_; ++s)
	{
		subindex_[s] = static_cast<int>(message[pos++]);
		subprimobj_[s] = message[pos++];
		subdualobj_[s] = message[pos++];
		CoinCopyN(message + pos, model_->getNumSubproblemCouplingCols(s), subsolution_[s]);
		pos += model_->getNumSubproblemCouplingCols(s);
	}
	worker_ = source;
	rcount_ = pos;

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
}

STO_RTN_CODE DdMaster::createMessage()
{
	BGN_TRY_CATCH

	int pos = 0;
	sbuf_[pos++] = static_cast<double>(status_);
	const double * theta = primsol_;
	const double * lambda = primsol_ + model_->getNumSubproblems();
	for (int s = 0; s < nsubprobs_; ++s)
	{
		int nlambda = model_->getNumSubproblemCouplingRows(s);
		const int * relevantRows = model_->getSubproblemCouplingRowIndices(s);
		/** theta */
		sbuf_[pos++] = theta[subindex_[s]];
		/** lambda */
		for (int k = 0; k < nlambda; ++k)
			sbuf_[pos++] = lambda[relevantRows[k]];
	}
	scount_ = pos;

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
}
