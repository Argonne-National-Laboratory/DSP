/*
 * DdMasterAtr.cpp
 *
 *  Created on: Mar 22, 2016
 *      Author: kibaekkim
 */

//#define DSP_DEBUG

#include "Solver/DualDecomp/DdMasterAtr.h"

DdMasterAtr::DdMasterAtr(
		DspParams *  par,
		DecModel *   model,
		DspMessage * message,
		int          nworkers,
		int          maxnumsubprobs):
	DdMasterTr(par, model, message, nworkers, maxnumsubprobs)
{
	/** number of cuts generated at the last iteration for each worker */
	nlastcuts_ = new int [nworkers_];
	CoinZeroN(nlastcuts_, nworkers_);

	/** solution for each worker process */
	primsol_to_worker_ = new double * [nworkers_];
}

DdMasterAtr::~DdMasterAtr()
{
	FREE_ARRAY_PTR(nlastcuts_);
	FREE_2D_ARRAY_PTR(nworkers_, primsol_to_worker_);
}

DSP_RTN_CODE DdMasterAtr::solve()
{
	BGN_TRY_CATCH

	DdMasterTr::solve();

	if (status_ == DSP_STAT_MW_CONTINUE)
	{
		for (unsigned i = 0; i < worker_.size(); ++i)
		{
			/** copy lambda */
			CoinCopyN(primsol_, nthetas_ + nlambdas_, primsol_to_worker_[worker_[i]]);
		}
//		DSPdebugMessage("primsol_:\n");
//		DSPdebug(StoMessage::printArray(nlambdas_+nthetas_, primsol_));
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMasterAtr::init()
{
	BGN_TRY_CATCH

	DdMasterTr::init();

	/** allocate memory for solution */
	for (int i = 0; i < nworkers_; ++i)
	{
		primsol_to_worker_[i] = new double [nthetas_ + nlambdas_];
		CoinFillN(primsol_to_worker_[i], nthetas_, COIN_DBL_MAX);
		CoinZeroN(primsol_to_worker_[i] + nthetas_, nlambdas_);
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMasterAtr::updateProblem()
{
	BGN_TRY_CATCH

	OsiCuts cuts;
	CoinPackedVector cutvec;
	double cutrhs = 0.0;

	for (unsigned i = 0; i < worker_.size(); ++i)
	{
//		DSPdebugMessage("worker %d\n", worker_[i]);
//		for (int j = 0, k = 0; j < nlambdas_ + nthetas_; ++j)
//		{
//			if (fabs(primsol_to_worker_[worker_[i]][j]) < 1.0e-10) continue;
//			if (k > 0 && k % 5 == 0) printf("\n");
//			printf("  [%6d] %+e", j, primsol_to_worker_[worker_[i]][j]);
//			k++;
//		}
//		printf("\n");
		/**
		 * Construct cuts of form
		 *   theta_s - (Hx - d) lambda_s <= D_s - (Hx - d) lambda_s_hat
		 */
		for (int s = 0; s < nsubprobs_[i]; ++s)
		{
//			StoMessage::printArray(model_->getNumSubproblemCouplingCols(subindex_[i][s]), subsolution_[i][s]);
			/** constructing */
			cutvec.clear();
			cutvec.insert(subindex_[i][s], 1.0); /**< theta part */
			cutrhs = subprimobj_[i][s];
			for (int j = 0; j < model_->getNumCouplingRows(); ++j)
			{
				/** evaluate solution on coupling constraints (if they are Hx = d, this is (Hx - d)_i) */
				double hx_d = model_->evalLhsCouplingRowSubprob(j, subindex_[i][s], subsolution_[i][s]) - model_->getRhsCouplingRow(j);
				if (fabs(hx_d) > 1.0e-10)
				{
					cutvec.insert(nthetas_ + j, -hx_d);
					cutrhs -= hx_d * primsol_to_worker_[worker_[i]][nthetas_ + j];
				}
			}

			/** cut placeholder */
			OsiRowCut * rc = new OsiRowCut;
			rc->setRow(cutvec);
			rc->setLb(-COIN_DBL_MAX);
			rc->setUb(cutrhs);
			rc->setEffectiveness(rc->violated(primsol_to_worker_[worker_[i]]));
			//DSPdebugMessage("cut violation %+e\n", rc->violated(primsol_to_worker_[worker_[i]]));

			if (rc->effectiveness() > 1.0e-6)
				/** local cut pool */
				cuts.insert(rc);
			else
				FREE_PTR(rc);
		}

		/** cut counter */
		nlastcuts_[worker_[i]] = cuts.sizeCuts();
	}
	message_->print(5, "-> master has %d rows, %d columns, and %d cuts to add\n",
			si_->getNumRows(), si_->getNumCols(), cuts.sizeCuts());

	/** add cut */
	if (cuts.sizeCuts() > 0)
	{
		//cuts.printCuts();
		si_->addCuts(cuts);

//		if (primobj_ < bestdualobj_)
//		{
//			stability_param_ *= 2;
//			message_->print(2, "-> update trust region size %+e\n", stability_param_);
//		}
	}
	else
	{
		/** number of cuts generated from the last iteration */
		int ncuts = 0;
		for (int i = 1; i < nworkers_; ++i)
			ncuts += nlastcuts_[i];

		if (ncuts == 0 && isSolutionBoundary() == true)
		{
			/** update trust region size */
			stability_param_ *= 2;
			message_->print(3, "-> update trust region size %+e\n", stability_param_);
		}
	}

	/** set trust region */
	if (status_ == DSP_STAT_MW_CONTINUE)
	{
		CoinCopyN(primsol_to_worker_[worker_[0]] + nthetas_, nlambdas_, stability_center_);
		setTrustRegion(stability_param_, stability_center_);
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMasterAtr::terminationTest()
{
	BGN_TRY_CATCH

	/** number of cuts generated from the last iteration */
	int ncuts = 0;
	for (int i = 1; i < nworkers_; ++i)
		ncuts += nlastcuts_[i];

	if (ncuts > 0)
		return status_;

	if (isSolutionBoundary() == false)
	{
		status_ = DSP_STAT_MW_STOP;
		message_->print(2, "-> STOP\n");
	}
	else
		message_->print(3, "The solution is at TR boundary.\n");

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return status_;
}
