/*
 * DdMasterReg.cpp
 *
 *  Created on: Feb 11, 2016
 *      Author: kibaekkim
 */

#include "Solver/DualDecomp/DdMasterReg.h"
#include "SolverInterface/SolverInterfaceOoqp.h"

DdMasterReg::DdMasterReg(DspParams * par, DecModel * model, StoMessage * message, int nworkers, int maxnumsubprobs):
	DdMaster(par, model, message, nworkers, maxnumsubprobs),
	nthetas_(0), nlambdas_(0),
	stability_param_(10.0), stability_center_(NULL), nlastcuts_(NULL),
	irowQ_(NULL), jcolQ_(NULL), dQ_(NULL), obj_(NULL) {}

DdMasterReg::~DdMasterReg()
{
	FREE_ARRAY_PTR(stability_center_);
	FREE_ARRAY_PTR(nlastcuts_);
	FREE_ARRAY_PTR(irowQ_);
	FREE_ARRAY_PTR(jcolQ_);
	FREE_ARRAY_PTR(dQ_);
	FREE_ARRAY_PTR(obj_);
}

/** initialize */
STO_RTN_CODE DdMasterReg::init()
{
	BGN_TRY_CATCH

	DdMaster::init();

	/** number of cuts generated at the last iteration for each worker */
	nlastcuts_ = new int [nworkers_];
	CoinZeroN(nlastcuts_, nworkers_);

	/** create problem */
	createProblem();

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
}

STO_RTN_CODE DdMasterReg::solve()
{
	BGN_TRY_CATCH

	if (status_ == STO_STAT_MW_STOP)
		return STO_RTN_OK;

	double cputime  = CoinCpuTime();
	double walltime = CoinGetTimeOfDay();

	/** solve */
	si_->solve();

	/** solver status */
	switch(si_->getStatus())
	{
	case STO_STAT_OPTIMAL:
	case STO_STAT_LIM_ITERorTIME:
	case STO_STAT_STOPPED_GAP:
	case STO_STAT_STOPPED_NODE:
	case STO_STAT_STOPPED_TIME:
	{
		/** get solution */
		CoinCopyN(si_->getSolution(), si_->getNumCols(), primsol_);
		//primobj_ = si_->getPrimalBound(); // This is not meaningful.
		primobj_ = 0.0;
		for (int j = 0; j < nthetas_; ++j)
			primobj_ += primsol_[j];

		/** update statistics */
		s_statuses_.push_back(si_->getStatus());
		s_primobjs_.push_back(si_->getPrimalBound());
		s_dualobjs_.push_back(si_->getDualBound());
		double * s_primsol = new double [si_->getNumCols()];
		CoinCopyN(si_->getSolution(), si_->getNumCols(), s_primsol);
		s_primsols_.push_back(s_primsol);
		s_primsol = NULL;
		s_cputimes_.push_back(CoinCpuTime() - cputime);
		s_walltimes_.push_back(CoinGetTimeOfDay() - walltime);

		break;
	}
	default:
		status_ = STO_STAT_MW_STOP;
		message_->print(0, "Warning: master solution status is %d\n", si_->getStatus());
		break;
	}

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
}

STO_RTN_CODE DdMasterReg::createProblem()
{
#define FREE_MEMORY       \
	FREE_ARRAY_PTR(ctype); \
	FREE_ARRAY_PTR(clbd); \
	FREE_ARRAY_PTR(cubd); \
	FREE_ARRAY_PTR(rlbd); \
	FREE_ARRAY_PTR(rubd); \
	FREE_ARRAY_PTR(bgn);  \
	FREE_ARRAY_PTR(len);  \
	FREE_ARRAY_PTR(ind);  \
	FREE_ARRAY_PTR(elem); \
	FREE_PTR(mat);

	int i, pos;
	int ncols, nrows, nzcnt;
	char * ctype = NULL;
	double * clbd = NULL;
	double * cubd = NULL;
	double * rlbd = NULL;
	double * rubd = NULL;
	CoinBigIndex * bgn = NULL;
	int * len = NULL;
	int * ind = NULL;
	double * elem = NULL;
	CoinPackedMatrix * mat = NULL;

	BGN_TRY_CATCH

	nthetas_  = model_->getNumSubproblems();
	nlambdas_ = model_->getNumCouplingRows();

	/** LP dimension */
	if (model_->nonanticipativity())
	{
		nrows = model_->getNumSubproblemCouplingCols(0); /** initial normalization constraint for nonanticipativity constraints */
		nzcnt = nrows * model_->getNumSubproblems();
	}
	else
	{
		nrows = 0;
		nzcnt = 0;
	}
	ncols = nthetas_ + nlambdas_;

	/** allocate memory */
	ctype = new char [ncols];
	clbd  = new double[ncols];
	cubd  = new double[ncols];
	rlbd  = new double[nrows];
	rubd  = new double[nrows];
	bgn   = new CoinBigIndex[nrows + 1];
	len   = new int[nrows];
	ind   = new int[nzcnt];
	elem  = new double[nzcnt];
	irowQ_ = new int [nlambdas_];
	jcolQ_ = new int [nlambdas_];
	dQ_    = new double [nlambdas_];
	obj_   = new double[ncols];
	stability_center_ = new double [nlambdas_];

	/** initial center */
	CoinZeroN(stability_center_, nlambdas_);

	/** all continuous variables */
	CoinFillN(ctype, ncols, 'C');

	/** c */
	CoinFillN(obj_, nthetas_, 1.0);
	CoinZeroN(obj_ + nthetas_, nlambdas_);

	/** bounds */
	CoinFillN(clbd, nthetas_, -COIN_DBL_MAX);
	CoinFillN(cubd, nthetas_, +COIN_DBL_MAX);

	/** nonnegative or nonpositive multipliers according to sense */
	for (i = 0; i < nlambdas_; i++)
	{
		if (model_->getSenseCouplingRow(i) == 'L')
			clbd[nthetas_ + i] = 0;
		else if (model_->getSenseCouplingRow(i) == 'G')
			cubd[nthetas_ + i] = 0;
	}

	if (model_->nonanticipativity())
	{
		/** row bounds */
		CoinZeroN(rlbd, nrows);
		CoinZeroN(rubd, nrows);

		/** for constraints */
		pos = 0;
		for (i = 0; i < nrows; ++i)
		{
			bgn[i] = pos;
			for (int j = 0; j < model_->getNumSubproblems(); ++j)
			{
				ind[pos] = nthetas_ + j * nrows + i;
				elem[pos] = 1.0;
				pos++;
			}
			len[i] = pos - bgn[i];
		}
		bgn[nrows] = pos;
		assert(pos == nzcnt);
	}

	/** constraint matrix */
	mat = new CoinPackedMatrix(false, ncols, nrows, nzcnt, elem, ind, bgn, len);
	DSPdebug(mat->verifyMtx(4));

	/** create solver interface */
	si_ = new SolverInterfaceOoqp(par_);

	/** [MAX]imization */
	si_->setObjSense(-1);

	/** copy problem data */
	si_->loadProblem(mat, clbd, cubd, obj_, ctype, rlbd, rubd);

	/** set lower triangle Hessian matrix */
	SolverInterfaceOoqp * ooqp = dynamic_cast<SolverInterfaceOoqp*>(si_);
	if (ooqp)
	{
		CoinIotaN(irowQ_, nlambdas_, nthetas_);
		CoinIotaN(jcolQ_, nlambdas_, nthetas_);
		CoinFillN(dQ_, nlambdas_, -0.5 / stability_param_);
		ooqp->setHessian(nlambdas_, irowQ_, jcolQ_, dQ_);
	}

	/** allocate memory for solution */
	primsol_ = new double [ncols];
	CoinZeroN(primsol_, ncols);

	/** set print level */
	si_->setPrintLevel(CoinMax(0, par_->getIntParam("LOG_LEVEL") - 2));

	END_TRY_CATCH_RTN(FREE_MEMORY,STO_RTN_ERR)

	FREE_MEMORY;

	return STO_RTN_OK;
#undef FREE_MEMORY
}

/** update problem */
STO_RTN_CODE DdMasterReg::updateProblem()
{
	OsiCuts cuts;
	CoinPackedVector cutvec;
	double cutrhs = 0.0;

	BGN_TRY_CATCH

	/**
	 * Construct cuts of form
	 *   theta_s - (Hx - d) lambda_s <= D_s - (Hx - d) lambda_s_hat
	 */
	for (int s = 0; s < nsubprobs_; ++s)
	{
		/** cut effectiveness */
		double eff = primsol_[subindex_[s]] - subprimobj_[s];
		if (eff < 1.0e-6) continue;

		/** constructing */
		cutvec.insert(subindex_[s], 1.0); /**< theta part */
		cutrhs = subprimobj_[s];
		for (int i = 0; i < model_->getNumCouplingRows(); ++i)
		{
			/** evaluate solution on coupling constraints (if they are Hx = d, this is (Hx - d)_i) */
			double hx_d = model_->evalLhsCouplingRowSubprob(i, subindex_[s], subsolution_[s]) - model_->getRhsCouplingRow(i);
			if (fabs(hx_d) > 1.0e-10)
			{
				cutvec.insert(nthetas_ + i, -hx_d);
				cutrhs -= hx_d * primsol_[nthetas_ + i];
			}
		}

		/** cut placeholder */
		OsiRowCut * rc = new OsiRowCut;
		rc->setRow(cutvec);
		rc->setLb(-COIN_DBL_MAX);
		rc->setUb(cutrhs);
		rc->setEffectiveness(eff);

		/** local cut pool */
		cuts.insert(rc);
	}

	/** cut counter */
	nlastcuts_[worker_] = cuts.sizeCuts();

	/** add cut */
	if (cuts.sizeCuts() > 0)
		si_->addCuts(cuts);
	else
	{
		/** number of cuts generated from the last iteration */
		int ncuts = 0;
		for (int i = 0; i < nworkers_; ++i)
			ncuts += nlastcuts_[i];

		if (ncuts == 0)
		{
			/** solve without regularization term */
			solveWithoutReg();

			/** update regularization term */
			stability_param_ *= 2;
		}
	}

	/** update objective function */
	for (int j = 0; j < nlambdas_; ++j)
	{
		obj_[nthetas_+j] = primsol_[nthetas_+j] / stability_param_;
		dQ_[j] = -0.5 / stability_param_;
	}
	refreshObjective();

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
}

/** update objective */
STO_RTN_CODE DdMasterReg::refreshObjective()
{
	BGN_TRY_CATCH

	SolverInterfaceOoqp * ooqp = dynamic_cast<SolverInterfaceOoqp*>(si_);
	ooqp->setObjCoef(obj_);
	ooqp->setHessian(nlambdas_, irowQ_, jcolQ_, dQ_);

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
}

/** solve without regularization */
STO_RTN_CODE DdMasterReg::solveWithoutReg()
{
	BGN_TRY_CATCH

	/** disable regularization term
	 * we don't need to backup the original objective,
	 * because it will be updated later in function updateProblem(); */
	CoinZeroN(obj_ + nthetas_, nlambdas_);
	CoinZeroN(dQ_, nlambdas_);
	refreshObjective();

	/** solve */
	si_->solve();

	/** solver status */
	switch(si_->getStatus())
	{
	case STO_STAT_OPTIMAL:
	case STO_STAT_LIM_ITERorTIME:
	case STO_STAT_STOPPED_GAP:
	case STO_STAT_STOPPED_NODE:
	case STO_STAT_STOPPED_TIME:
		if (si_->getPrimalBound() - primobj_ < 1.0e-8)
		{
			/** let's stop */
			status_ = STO_STAT_MW_STOP;
		}
		break;
	default:
		message_->print(0, "Warning: master solution status is %d\n", si_->getStatus());
		break;
	}

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
}
