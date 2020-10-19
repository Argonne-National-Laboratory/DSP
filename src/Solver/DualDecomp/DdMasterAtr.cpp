/*
 * DdMasterAtr.cpp
 *
 *  Created on: Mar 22, 2016
 *      Author: kibaekkim
 */

//#define DSP_DEBUG

#include "Solver/DualDecomp/DdMasterAtr.h"
#include "SolverInterface/DspOsiOoqpEps.h"

DdMasterAtr::DdMasterAtr(
		DecModel *   model,   /**< model pointer */
		DspParams *  par,     /**< parameter pointer */
		DspMessage * message, /**< message pointer */
		int nworkers          /**< number of workers */):
DdMasterTr(model, par, message),
nworkers_(nworkers),
is_updated_(false)
{
	/** number of cuts generated at the last iteration for each worker */
	nlastcuts_ = new int [nworkers_];
	CoinZeroN(nlastcuts_, nworkers_);

	/** solution for each worker process */
	primsol_to_worker_ = new double * [nworkers_];

	proved_optimality_ = new bool [nworkers_];
}

DdMasterAtr::DdMasterAtr(const DdMasterAtr& rhs) :
DdMasterTr(rhs),
nworkers_(rhs.nworkers_),
worker_(rhs.worker_),
solution_key_(rhs.solution_key_),
nsubprobs_(rhs.nsubprobs_),
subindex_(rhs.subindex_),
subprimobj_(rhs.subprimobj_),
subdualobj_(rhs.subdualobj_),
subsolution_(rhs.subsolution_),
is_updated_(rhs.is_updated_) {
	/** number of cuts generated at the last iteration for each worker */
	nlastcuts_ = new int [nworkers_];
	CoinCopyN(rhs.nlastcuts_, nworkers_, nlastcuts_);

	/** solution for each worker process */
	primsol_to_worker_ = new double * [nworkers_];
	for (int i = 0; i < nworkers_; ++i) {
		primsol_to_worker_[i] = new double [nthetas_ + nlambdas_];
		CoinCopyN(rhs.primsol_to_worker_[i], nthetas_+nlambdas_, primsol_to_worker_[i]);
	}

	/** indication of the proof of optimality */
	proved_optimality_ = new bool [nworkers_];
	CoinCopyN(rhs.proved_optimality_, nworkers_, proved_optimality_);
}

DdMasterAtr::~DdMasterAtr()
{
	FREE_ARRAY_PTR(nlastcuts_);
	FREE_ARRAY_PTR(proved_optimality_);
	FREE_2D_ARRAY_PTR(nworkers_, primsol_to_worker_);
}

DSP_RTN_CODE DdMasterAtr::solve()
{
	BGN_TRY_CATCH

	DSP_RTN_CHECK_RTN_CODE(DdMasterTr::solve());

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMasterAtr::init()
{
	BGN_TRY_CATCH

	DdMasterTr::init();

	/** initial primal solution */
	for (int i = 0; i < nworkers_; ++i)
	{
		primsol_to_worker_[i] = new double [nthetas_ + nlambdas_];
		CoinFillN(primsol_to_worker_[i], nthetas_, COIN_DBL_MAX);
		CoinZeroN(primsol_to_worker_[i] + nthetas_, nlambdas_);
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMasterAtr::updateProblem(
		double * primsol, /**< master primal solution at which newdual is obtained */
		double   newdual  /**< new dual objective value */)
{
	BGN_TRY_CATCH
#if 0
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
		double * primsol = master_primsol_[i];
		for (int s = 0; s < nsubprobs_[i]; ++s)
		{
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
					cutrhs -= hx_d * primsol[nthetas_ + j];
				}
			}

			/** cut placeholder */
			OsiRowCut * rc = new OsiRowCut;
			rc->setRow(cutvec);
			rc->setLb(-COIN_DBL_MAX);
			rc->setUb(cutrhs);
			rc->setEffectiveness(rc->violated(primsol));
			//DSPdebugMessage("cut violation %+e\n", rc->violated(primsol));

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
			getSiPtr()->getNumRows(), getSiPtr()->getNumCols(), cuts.sizeCuts());

	/** add cut */
	if (cuts.sizeCuts() > 0)
	{
		//cuts.printCuts();
		getSiPtr()->addCuts(cuts);
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
			message_->print(1, "-> update trust region size %+e\n", stability_param_);
		}
	}

	/** set trust region */
	if (status_ == DSP_STAT_MW_CONTINUE)
	{
		CoinCopyN(primsol + nthetas_, nlambdas_, stability_center_);
		setTrustRegion(stability_param_, stability_center_);
	}
#endif

	/** retrieve master primal solution */
	int nCutsAdded = 0;
	double curprimobj = 0.0;//getSiPtr()->getPrimalBound(); /** current primal objective value */
	for (int s = 0; s < nthetas_; ++s)
		curprimobj += primsol[s];

	/** TODO calculate primal/dual objectives */
#if 0
	double newprimal = 0.0;
	double newdual = 0.0;
	for (int s = 0; s < model_->getNumSubproblems(); ++s)
	{
		newprimal += subprimobj_[s];
		newdual += subdualobj_[s];
	}
#endif

	if (newdual >= bestdualobj_ + 1.0e-4 * (curprimobj - bestdualobj_))
	{
		message_->print(2, "TR  %s STEP: dual objective %e", isSolved_ ? "SERIOUS" : "INITIAL", newdual);

		/** mark cuts not to be deleted */
		for (int i = cuts_->sizeCuts() - nCutsAdded; i < cuts_->sizeCuts(); ++i)
			possiblyDelete_[i] = false;

		if (isSolved_)
		{
			/** update proximal point */
			CoinCopyN(primsol + nthetas_, nlambdas_, stability_center_);
			message_->print(3, ", updated proximal point");

			/** possibly delete cuts */
			//possiblyDeleteCuts(newdual);

			/** is solution boundary? */
			if (isSolutionBoundary() &&
				/*primalBound - dual_bound < 0 ||*/
				newdual >= bestdualobj_ + 0.5 * (curprimobj - bestdualobj_))
			{
				/** increase trust region */
				stability_param_ = CoinMin(2. * stability_param_, 1.0e+4);
				message_->print(3, ", increased trust region size %e", stability_param_);
			}

			/** set trust region */
			setTrustRegion(stability_param_, stability_center_);
		}

		/** update dual bound */
		bestdualobj_ = newdual;
		trcnt_ = 0;

		message_->print(2, "\n");
	}
	else
	{
		/** add cuts and increase minor cut counter */
		ncuts_minor_ += addCuts();

		/** null step */
		message_->print(3, "[TR]  NULL STEP: dual objective %e", newdual);

		/** The following rule is from Linderoth and Wright (2003) */
		int nullsteps_allowed = 3;
		double rho = CoinMin(1.0, stability_param_) * CoinMax(bestdualobj_ - newdual, linerr_) / (curprimobj - bestdualobj_);
		message_->print(3, ", rho %e, trcnt %d", rho, trcnt_);
		if (rho > 0) trcnt_++;
		if (rho >= 3 || (trcnt_ >= nullsteps_allowed && fabs(rho - 2.) < 1.0))
		{
			/** decrease trust region */
			stability_param_ = CoinMax(1.0e-2, stability_param_/CoinMin(rho, 4.));
			message_->print(3, ", decreased trust region size %e", stability_param_);
			trcnt_ = 0;

			/** set trust region */
			setTrustRegion(stability_param_, stability_center_);
		}
		message_->print(3, "\n");
	}

#ifdef DSP_HAS_OOQP
	DspOsiOoqpEps * ooqp = dynamic_cast<DspOsiOoqpEps*>(osi_);
	if (ooqp)
	{
		if (ooqp->ooqp_->hasOoqpStatus_ && isSolved_)
		{
			DSPdebugMessage("bestprimobj %+e bestdualobj %+e\n", bestprimobj_, bestdualobj_);
			double epsilon = (getSiPtr()->getObjValue() - newdual + ooqp->ooqp_->getDualityGap()) / (1. + fabs(getSiPtr()->getObjValue()));
			if (epsilon > 1.) epsilon = 1.;
			ooqp->ooqp_->setOoqpStatus(epsilon, -bestprimobj_, -bestdualobj_);
		}
	}
#endif

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMasterAtr::updateTrustRegion(const double * primsol)
{
	BGN_TRY_CATCH

	/** add cuts and increase minor cut counter */
	ncuts_minor_ += addCuts();
#if 0
	/** current primal objective value */
	double curprimobj = 0.0;
	for (int s = 0; s < nthetas_; ++s)
		curprimobj += primsol[s];

	/** The following rule is from Linderoth and Wright (2003) */
	int nullsteps_allowed = 3 * par_->getIntPtrParamSize("ARR_PROC_IDX");
	double rho = CoinMin(1.0, stability_param_) * linerr_ / (curprimobj - bestdualobj_);
	message_->print(3, ", rho %e, trcnt %d", rho, trcnt_);
	if (rho > 0) trcnt_++;
	if (rho >= 3 || (trcnt_ >= nullsteps_allowed && fabs(rho - 2.) < 1.0))
	{
		/** decrease trust region */
		stability_param_ *= 1.0 / CoinMin(rho, 4.);
		message_->print(3, ", decreased trust region size %e", stability_param_);
		trcnt_ = 0;

		/** set trust region */
		setTrustRegion(stability_param_, stability_center_);
	}
	message_->print(3, "\n");
#endif
	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMasterAtr::terminationTest()
{
	int signal = status_;

	BGN_TRY_CATCH

#if 0
	/** number of cuts generated from the last iteration */
	int ncuts = 0;
	for (int i = 1; i < nworkers_; ++i)
		ncuts += nlastcuts_[i];

	if (ncuts > 0)
		return status_;
#endif

	if (status_ == DSP_STAT_MW_STOP) return status_;

	double absgap = primobj_ - bestdualobj_;
	double relgap = fabs(absgap) / (1.e-10 + fabs(primobj_));
	if (relgap <= par_->getDblParam("DD/STOP_TOL") + par_->getDblParam("DD/SUB/GAPTOL"))
	{
		signal = DSP_STAT_MW_STOP;
		status_ = DSP_STAT_OPTIMAL;
		message_->print(0, "TR  STOP with gap tolerance %+e (%.2f%%).\n", absgap, relgap*100);
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return signal;
}

int DdMasterAtr::addCuts(bool possiblyDel)
{
	OsiCuts cuts;
	CoinPackedVector cutvec;
	double cutrhs = 0.0;
	int nCutsAdded = 0;

	BGN_TRY_CATCH

	/** initialize linearization error */
	int ncuts = 0;
	linerr_ = 0.0;

//	DSPdebugMessage("Number of workers %lu\n", worker_.size());
	for (unsigned i = 0; i < worker_.size(); ++i)
	{
		/** retrieve master primal solution */
		double * primsol = primsol_to_worker_[worker_[i]-1];//master_primsol_[i];
		double primbound = 0.0;
		for (int s = 0; s < nthetas_; ++s)
			primbound += primsol[s];

		/**
		 * Construct cuts of form
		 *   theta_s - (Hx - d) lambda_s <= D_s - (Hx - d) lambda_s_hat
		 */
		for (int s = 0; s < nsubprobs_[i]; ++s)
		{
			/** constructing cut and calculating error*/
			cutvec.clear();
			cutvec.insert(subindex_[i][s], 1.0); /**< theta part */
			linerr_ += subprimobj_[i][s];
			cutrhs = subprimobj_[i][s];
			for (int j = 0; j < model_->getNumCouplingRows(); ++j)
			{
				/** evaluate solution on coupling constraints (if they are Hx = d, this is (Hx - d)_i) */
				double hx_d = model_->evalLhsCouplingRowSubprob(j, subindex_[i][s], subsolution_[i][s]) - model_->getRhsCouplingRow(j);
				if (fabs(hx_d) > 1.0e-10)
				{
					cutvec.insert(nthetas_ + j, -hx_d);
					cutrhs -= hx_d * primsol[nthetas_ + j];
					linerr_ += hx_d * (stability_center_[j] - primsol[nthetas_ + j]);
				}
			}

			/** count number of cuts generated */
			ncuts++;

			/** cut placeholder */
			OsiRowCut * rc = new OsiRowCut;
			rc->setRow(cutvec);
			rc->setLb(-COIN_DBL_MAX);
			rc->setUb(cutrhs);
			rc->setEffectiveness(rc->violated(primsol));
			//message_->print(2, "cut violation %+e\n", rc->effectiveness());
			//DSPdebugMessage("cut violation %+e\n", rc->violated(primsol));

			if (rc->effectiveness() > 1.0e-6)
			{
//#define ADD_ALL_CUTS
#ifndef ADD_ALL_CUTS
				/** number of cuts before adding cut */
				int nCutsBefore = cuts_->sizeCuts();

				/** add cut if not duplicate */
				cuts_->insertIfNotDuplicate(*rc);

				if (nCutsBefore < cuts_->sizeCuts())
				{
					/** insertIfNotDuplicate does not set effectiveness */
					cuts_->rowCutPtr(nCutsBefore)->setEffectiveness(rc->effectiveness());
					cuts_age_.push_back(0);
					possiblyDelete_.push_back(possiblyDel);
					masterobjsAtCutAdd_.push_back(primbound);
#endif
					/** local cut pool */
					cuts.insert(rc);
#ifndef ADD_ALL_CUTS
				}
#endif
			}
			else
				FREE_PTR(rc);
		}

		/** cut counter */
		nlastcuts_[worker_[i]-1] = cuts.sizeCuts();
	}
	message_->print(5, "-> master has %d rows, %d columns, and %d cuts to add\n",
			getSiPtr()->getNumRows(), getSiPtr()->getNumCols(), cuts.sizeCuts());

	/** calculate linearization error */
	linerr_ -= bestdualobj_;
	linerr_ *= model_->getNumSubproblems() * 1.0 / ncuts;

	/** add cut */
	nCutsAdded = cuts.sizeCuts();
	if (nCutsAdded > 0)
	{
		//cuts.printCuts();
		getSiPtr()->applyCuts(cuts);
		is_updated_ = true;
		CoinFillN(proved_optimality_, nworkers_, false);
	}
	else
	{
		/** TODO recruit back some cuts if no cut is generated */
		message_->print(10, "WARNING: Cut recruit is not implemented.\n");
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return nCutsAdded;
}

DSP_RTN_CODE DdMasterAtr::clearSubprobData()
{
	BGN_TRY_CATCH

	worker_.clear();
	solution_key_.clear();
	nsubprobs_.clear();
	subindex_.clear();
	subprimobj_.clear();
	subdualobj_.clear();
	subsolution_.clear();

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMasterAtr::setPrimsolToWorker(
		int      worker_id, /**< worker ID */
		double * primsol    /**< primal solution assigned to worker */)
{
	BGN_TRY_CATCH

	CoinCopyN(primsol, getSiPtr()->getNumCols(), primsol_to_worker_[worker_id]);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMasterAtr::setTrustRegion(
		double stability_param,
		double* stability_center)
{
	BGN_TRY_CATCH

	DdMasterTr::setTrustRegion(stability_param, stability_center);
	is_updated_ = true;
	CoinFillN(proved_optimality_, nworkers_, false);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}
