/*
 * DdWorkerUB.cpp
 *
 *  Created on: Mar 28, 2016
 *      Author: kibaekkim
 */

//#define DSP_DEBUG
#include <Solver/DualDecomp/DdWorkerUB.h>

DdWorkerUB::DdWorkerUB(DspParams * par, DecModel * model, DspMessage * message):
DdWorkerLB(par, model, message)
{
	// TODO Auto-generated constructor stub

}

DdWorkerUB::~DdWorkerUB() {
	// TODO Auto-generated destructor stub
}

DSP_RTN_CODE DdWorkerUB::fixCouplingVariableValues(CoinPackedVector * val)
{
	BGN_TRY_CATCH

	for (unsigned s = 0; s < subprobs_.size(); ++s)
	{
		int ncols = model_->getNumSubproblemCouplingCols(subprobs_[s]->sind_);
		double * denseval = val->denseVector(ncols);
		for (int j = 0; j < ncols; ++j)
			subprobs_[s]->si_->setColBounds(j, denseval[j], denseval[j]);
		FREE_ARRAY_PTR(denseval);
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdWorkerUB::solve() {
	double cputime;
	double walltime;

	BGN_TRY_CATCH

	double primobj = 0.0;
	double dualobj = 0.0;
	double total_cputime = 0.0;
	double total_walltime = 0.0;
	double pargaptol = par_->getDblParam("DD/STOP_TOL");

	for (unsigned s = 0; s < subprobs_.size(); ++s) {

		cputime = CoinCpuTime();
		walltime = CoinGetTimeOfDay();

		/** update time stamp */
		ticToc();

		/** set time limit */
		subprobs_[s]->setTimeLimit(
				CoinMin(time_remains_,
						par_->getDblParam("SCIP/TIME_LIM")));
		/** solve */
		subprobs_[s]->solve();

		/** check status. there might be unexpected results. */
		switch (subprobs_[s]->si_->getStatus()) {
		case DSP_STAT_OPTIMAL:
		case DSP_STAT_LIM_ITERorTIME:
		case DSP_STAT_STOPPED_GAP:
		case DSP_STAT_STOPPED_NODE:
		case DSP_STAT_STOPPED_TIME:
			break;
		default:
			status_ = DSP_STAT_MW_STOP;
			message_->print(0,
					"Warning: subproblem %d solution status is %d\n", s,
					subprobs_[s]->si_->getStatus());
			break;
		}
		if (status_ == DSP_STAT_MW_STOP)
			break;

		primobj += subprobs_[s]->si_->getPrimalBound();
		dualobj += subprobs_[s]->si_->getDualBound();
		total_cputime += CoinCpuTime() - cputime;
		total_walltime += CoinGetTimeOfDay() - walltime;
	}

	/** update statistics */
	s_statuses_.push_back(status_);
	s_primobjs_.push_back(primobj);
	s_dualobjs_.push_back(dualobj);
	s_cputimes_.push_back(total_cputime);
	s_walltimes_.push_back(total_walltime);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}
