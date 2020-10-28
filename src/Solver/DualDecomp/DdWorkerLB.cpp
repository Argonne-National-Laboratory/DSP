/*
 * DdWorkerLB.cpp
 *
 *  Created on: Apr 4, 2016
 *      Author: kibaekkim
 */

//#define DSP_DEBUG

#include "Solver/DualDecomp/DdWorkerLB.h"

DdWorkerLB::DdWorkerLB(
		DecModel *   model,  /**< model pointer */
		DspParams *  par,    /**< parameter pointer */
		DspMessage * message /**< message pointer */) :
DdWorker(model, par, message), 
solution_key_(-1), 
isInit_(true) {}

DdWorkerLB::DdWorkerLB(const DdWorkerLB& rhs) :
DdWorker(rhs), 
solution_key_(rhs.solution_key_), 
isInit_(rhs.isInit_) {
	for (unsigned s = 0; s < rhs.subprobs_.size(); ++s)
		subprobs_[s] = rhs.subprobs_[s]->clone();
}

DdWorkerLB::~DdWorkerLB() {
	for (unsigned s = 0; s < subprobs_.size(); ++s)
		FREE_PTR(subprobs_[s]);
}

DSP_RTN_CODE DdWorkerLB::init() {
	BGN_TRY_CATCH
	/** status */
	status_ = DSP_STAT_MW_CONTINUE;
	/** create problem */
	DSP_RTN_CHECK_THROW(createProblem(
		par_->getIntPtrParamSize("ARR_PROC_IDX"), par_->getIntPtrParam("ARR_PROC_IDX")));
	END_TRY_CATCH_RTN(;, DSP_RTN_ERR)
	return DSP_RTN_OK;
}

DSP_RTN_CODE DdWorkerLB::solve() {
    double cputime;
    double walltime;
	int rtncode = DSP_RTN_OK;

    BGN_TRY_CATCH

	double primobj = 0.0;
	double dualobj = 0.0;
	double total_cputime = 0.0;
	double total_walltime = 0.0;
	double pargaptol = par_->getDblParam("DD/STOP_TOL");

	for (unsigned s = 0; s < subprobs_.size(); ++s) {
		cputime = CoinCpuTime();
		walltime = CoinGetTimeOfDay();

		/** reset gap tolerance */
		if (subprobs_[s]->getSiPtr()->getNumIntegers() > 0)
			subprobs_[s]->setGapTol(par_->getDblParam("DD/SUB/GAPTOL"));

		bool resolve = true;
		while (resolve) {
			resolve = false;

			/** set time limit */
			subprobs_[s]->setTimeLimit(
					CoinMin(CoinMax(0.01, time_remains_), par_->getDblParam("DD/SUB/TIME_LIM")));
			/** solve */
			subprobs_[s]->solve();

			/** check status. there might be unexpected results. */
			switch (subprobs_[s]->getStatus()) {
				case DSP_STAT_STOPPED_TIME:
				case DSP_STAT_LIM_ITERorTIME:
				case DSP_STAT_LIM_INFEAS:
				case DSP_STAT_STOPPED_GAP:
				case DSP_STAT_STOPPED_NODE:
					message_->print(3, "Warning: subproblem %d solution status is %d\n",
									subprobs_[s]->sind_, subprobs_[s]->getStatus());
					break;
				case DSP_STAT_OPTIMAL:
					break;
				default:
					status_ = DSP_STAT_MW_STOP;
					rtncode = DSP_RTN_ERR;
					message_->print(0, "Warning: subproblem %d solution status is %d\n",
									subprobs_[s]->sind_, subprobs_[s]->getStatus());
					break;
			}
			if (status_ == DSP_STAT_MW_STOP)
				break;

			/** set solution gap tolerance */
			if (subprobs_[s]->getSiPtr()->getNumIntegers() > 0) {
				if (isInit_ == false &&
						subprobs_[s]->getPrimalObjective() >= subprobs_[s]->theta_ &&
						subprobs_[s]->gapTol_ > pargaptol) {
					/** TODO parameterize this */
					double gapTol = subprobs_[s]->gapTol_ * 0.5;
					if (gapTol < pargaptol)
						gapTol = pargaptol;
					subprobs_[s]->setGapTol(gapTol);
					resolve = true;
				}
			}
		}

		message_->print(10, "-> subprob %d: solved with gap tolerance %e \n",
						subprobs_[s]->sind_, subprobs_[s]->gapTol_);

		if (status_ == DSP_STAT_MW_STOP)
			break;

		if (subprobs_[s]->getStatus() == DSP_STAT_LIM_INFEAS)
			primobj = COIN_DBL_MAX;
		else if (primobj < COIN_DBL_MAX)
			primobj += subprobs_[s]->getDspOsiPtr()->getPrimObjValue();
		dualobj += subprobs_[s]->getDspOsiPtr()->getDualObjValue();
		total_cputime += CoinCpuTime() - cputime;
		total_walltime += CoinGetTimeOfDay() - walltime;

		/** consume time */
		time_remains_ -= CoinGetTimeOfDay() - walltime;
	}

	isInit_ = false;

	/** update statistics */
	s_statuses_.push_back(status_);
	s_primobjs_.push_back(primobj);
	s_dualobjs_.push_back(dualobj);
	s_cputimes_.push_back(total_cputime);
	s_walltimes_.push_back(total_walltime);

    END_TRY_CATCH_RTN(;, DSP_RTN_ERR)

    return rtncode;
}

DSP_RTN_CODE DdWorkerLB::createProblem(int nsubprobs, int* subindex) {
	BGN_TRY_CATCH

	/** release all the subproblem data */
	for (unsigned s = 0; s < subprobs_.size(); ++s)
		FREE_PTR(subprobs_[s]);
	subprobs_.clear();

	for (int s = 0; s < nsubprobs; ++s) {
        /** create subproblem instance */
        DdSub *subprob = new DdSub(subindex[s], par_, model_, message_);
        /** initialize */
        DSP_RTN_CHECK_THROW(subprob->init());
        assert(subprob->getSiPtr());
        /** store */
        subprobs_.push_back(subprob);
    }
	END_TRY_CATCH_RTN(;, DSP_RTN_ERR)
	return DSP_RTN_OK;
}

