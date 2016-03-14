/*
 * DdWorker.cpp
 *
 *  Created on: Feb 9, 2016
 *      Author: kibaekkim
 */

/** DSP */
#include "Solver/DualDecomp/DdWorker.h"
#include "SolverInterface/SolverInterfaceScip.h"

DdWorker::DdWorker(DspParams * par, DecModel * model, StoMessage * message) :
	DdSolver(par, model, message) {}

DdWorker::~DdWorker()
{
	for (unsigned s = 0; s < subprobs_.size(); ++s)
		FREE_PTR(subprobs_[s]);
}

STO_RTN_CODE DdWorker::init()
{
	BGN_TRY_CATCH

	/** message count */
	scount_ = 1;
	rcount_ = 1;
	for (unsigned s = 0; s < par_->getIntPtrParamSize("ARR_PROC_IDX"); ++s)
	{
		int ss = par_->getIntPtrParam("ARR_PROC_IDX")[s];
		scount_ += 3 + model_->getNumSubproblemCouplingCols(ss);
		rcount_ += 1 + model_->getNumSubproblemCouplingRows(ss);
	}
	/** message buffer */
	sbuf_ = new double [scount_];
	sbuf_[0] = static_cast<double>(par_->getIntPtrParamSize("ARR_PROC_IDX"));

	/** status */
	status_ = STO_STAT_MW_CONTINUE;

	/** create problem */
	createProblem();

	/** time stamp */
	ticToc();

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
}

STO_RTN_CODE DdWorker::solve()
{
	double cputime;
	double walltime;

	BGN_TRY_CATCH

	double primobj = 0.0;
	double dualobj = 0.0;
	double total_cputime  = 0.0;
	double total_walltime = 0.0;

	for (unsigned s = 0; s < subprobs_.size(); ++s)
	{
		cputime  = CoinCpuTime();
		walltime = CoinGetTimeOfDay();

		bool resolve = true;
		while(resolve)
		{
			resolve = false;

			/** update time stamp */
			ticToc();
			/** set time limit */
			subprobs_[s]->setTimeLimit(CoinMin(time_remains_, par_->getDblParam("SCIP/TIME_LIM")));
			/** solve */
			subprobs_[s]->solve();

			/** check status. there might be unexpected results. */
			switch(subprobs_[s]->si_->getStatus())
			{
			case STO_STAT_OPTIMAL:
			case STO_STAT_LIM_ITERorTIME:
			case STO_STAT_STOPPED_GAP:
			case STO_STAT_STOPPED_NODE:
			case STO_STAT_STOPPED_TIME:
				break;
			default:
				status_ = STO_STAT_MW_STOP;
				message_->print(0, "Warning: subproblem %d solution status is %d\n", s, subprobs_[s]->si_->getStatus());
				break;
			}
			if (status_ == STO_STAT_MW_STOP)
				break;

			/** set solution gap tolerance */
			if (subprobs_[s]->getPrimalBound() >= subprobs_[s]->theta_
					&& subprobs_[s]->gapTol_ > par_->getDblParam("SCIP/GAP_TOL"))
			{
				/** TODO parameterize this */
				double gapTol = subprobs_[s]->gapTol_ * 0.5;
				if (gapTol < par_->getDblParam("SCIP/GAP_TOL"))
					gapTol = par_->getDblParam("SCIP/GAP_TOL");
				subprobs_[s]->setGapTol(gapTol);
				resolve = true;
			}
		}

		message_->print(2, "-> subprob %d: solved with gap tolerance %e \n", subprobs_[s]->sind_, subprobs_[s]->gapTol_);

		if (status_ == STO_STAT_MW_STOP)
			break;

		primobj += subprobs_[s]->si_->getPrimalBound();
		dualobj += subprobs_[s]->si_->getDualBound();
		total_cputime  += CoinCpuTime() - cputime;
		total_walltime += CoinGetTimeOfDay() - walltime;

		/** update statistics */
//		subprobs_[s]->s_statuses_.push_back(subprobs_[s]->si_->getStatus());
//		subprobs_[s]->s_primobjs_.push_back(subprobs_[s]->si_->getPrimalBound());
//		subprobs_[s]->s_dualobjs_.push_back(subprobs_[s]->si_->getDualBound());
//		double * s_primsol = new double [subprobs_[s]->si_->getNumCols()];
//		CoinCopyN(subprobs_[s]->si_->getSolution(), subprobs_[s]->si_->getNumCols(), s_primsol);
//		subprobs_[s]->s_primsols_.push_back(s_primsol);
//		s_primsol = NULL;
//		subprobs_[s]->s_cputimes_.push_back(CoinCpuTime() - cputime);
//		subprobs_[s]->s_walltimes_.push_back(CoinGetTimeOfDay() - walltime);
	}

	/** update statistics */
	s_statuses_.push_back(status_);
	s_primobjs_.push_back(primobj);
	s_dualobjs_.push_back(dualobj);
	s_cputimes_.push_back(total_cputime);
	s_walltimes_.push_back(total_walltime);

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
}

STO_RTN_CODE DdWorker::recvMessage(int size, double* message)
{
	if (size != rcount_)
	{
		message_->print(0, "Received message count is invalid.\n");
		return STO_RTN_ERR;
	}

	BGN_TRY_CATCH

	/** always message[0] == STO_STAT_MW_CONTINUE */
	assert(message[0] == STO_STAT_MW_CONTINUE);
	/** parse message */
	int pos = 1;
	for (unsigned s = 0; s < subprobs_.size(); ++s)
	{
		subprobs_[s]->theta_ = message[pos++];
		subprobs_[s]->updateProblem(message + pos);
		pos += model_->getNumSubproblemCouplingRows(s);
	}
	rcount_ = pos;

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
}

STO_RTN_CODE DdWorker::createProblem()
{
	BGN_TRY_CATCH

	for (int s = 0; s < par_->getIntPtrParamSize("ARR_PROC_IDX"); ++s)
	{
		/** create subproblem instance */
		message_->print(999, "Creating DD subproblem for index %d\n", par_->getIntPtrParam("ARR_PROC_IDX")[s]);
		DdSub * subprob = new DdSub(par_->getIntPtrParam("ARR_PROC_IDX")[s], par_);
		subprobs_.push_back(subprob);

		/** create subproblem */
		subprob->createProblem(model_);
		assert(subprob->si_);

		/** TODO general solver */
		/** set display level */
		SolverInterfaceScip * si = dynamic_cast<SolverInterfaceScip*>(subprob->si_);
		int loglevel = CoinMax(0, par_->getIntParam("LOG_LEVEL") - 3);
		if (si)
		 	subprob->si_->setPrintLevel(CoinMin(loglevel, 5));
		else
		 	subprob->si_->setPrintLevel(loglevel);
	}

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
}

/** create message */
STO_RTN_CODE DdWorker::createMessage()
{
	BGN_TRY_CATCH

	int pos = 1;
	for (unsigned s = 0; s < subprobs_.size(); ++s)
	{
		sbuf_[pos++] = subprobs_[s]->sind_;
		sbuf_[pos++] = subprobs_[s]->getPrimalBound();
		sbuf_[pos++] = subprobs_[s]->getDualBound();
		CoinCopyN(subprobs_[s]->si_->getSolution(), subprobs_[s]->ncols_coupling_, sbuf_ + pos);
		pos += subprobs_[s]->ncols_coupling_;
	}
	assert(pos == scount_);

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
}
