/*
 * DwBundleDualPips.cpp
 *
 *  Created on: Dec 7, 2017
 *      Author: kibaekkim
 */

#define DSP_DEBUG

#include "Solver/DantzigWolfe/DwBundleDualPips.h"
#include "Solver/DantzigWolfe/DwWorkerPips.h"
#include "Utility/DspUtility.h"

#define DEV071118

DwBundleDualPips& DwBundleDualPips::operator =(const DwBundleDualPips& rhs) {
	DwBundleDualSmip::operator =(rhs);
	pips_worker_ = rhs.pips_worker_;
	return *this;
}

void DwBundleDualPips::initDualSolver(
		const CoinPackedMatrix& m, 
		std::vector<double>& clbd, 
		std::vector<double>& cubd, 
		std::vector<double>& obj, 
		std::vector<double>& rlbd, 
		std::vector<double>& rubd) {
	/** does nothing */
}

DSP_RTN_CODE DwBundleDualPips::updateCenter(double penalty) {
	u_ = penalty;
	return DSP_RTN_OK;
}

DSP_RTN_CODE DwBundleDualPips::callMasterSolver() {
	BGN_TRY_CATCH

	/** sync columns (bundle information) only */
	DSP_RTN_CHECK_RTN_CODE(pips_worker_->sync(bestdualsol_, cols_generated_));

	/** run PIPS-IPM */
	DSP_RTN_CHECK_RTN_CODE(pips_worker_->solvePips(u_));

	/** TODO: parse PIPS-IPM solution status */
	status_ = DSP_STAT_OPTIMAL;

	END_TRY_CATCH_RTN(;,DSP_STAT_ABORT)

	return status_;
}

void DwBundleDualPips::assignMasterSolution(std::vector<double>& sol) {
	sol = bestdualsol_;
	for (int s = 0; s < tss_->getNumScenarios(); ++s) {
		sol[s] = pips_worker_->pips_theta_[s][0];
		assert(pips_worker_->pips_beta_[s].size() == tss_->getNumCols(0));
		for (size_t j = 0; j < tss_->getNumCols(0); ++j)
			sol[tss_->getNumScenarios() + tss_->getNumCols(0) * s + j] += pips_worker_->pips_beta_[s][j];
	}
}

double DwBundleDualPips::getObjValue() {
	return pips_worker_->pips_objval_ * -1.0;
}

void DwBundleDualPips::addDualRow(const CoinPackedVector& v, double lb, double ub) {	
	/** does nothing */
}

void DwBundleDualPips::removeAllDualRows() {	
	/** does nothing */
}

void DwBundleDualPips::printIterInfo() {
	message_->print(1, "Iteration %3d: DW Bound %+e, ", itercnt_, primobj_);
	message_->print(3, "Dual %+e, Approx %+e, ", -dualobj_, -bestdualobj_-v_);
	message_->print(1, "Best Dual %+e ", -bestdualobj_);
	if (relgap_ < 1000)
		message_->print(1, "(gap %.2f %%), ", relgap_*100);
	else
		message_->print(1, "(gap Large %%), ");
	// message_->print(1, "nrows %d, ncols %d, ", pips->getNumRows(), pips->getNumCols());
	message_->print(1, "timing (total %.2f, master %.2f, gencols %.2f), statue %d\n",
			t_total_, t_master_, t_colgen_, status_);
	message_->print(3, "  predicted ascent %+e, |p| %+e, alpha %+e, linerr %+e, eps %+e, u %+e, counter %d\n",
			-v_, absp_, alpha_, -linerr_, eps_, u_, counter_);

	/** log */
	log_time_.push_back(CoinGetTimeOfDay());
	log_bestdual_bounds_.push_back(-bestdualobj_);
}
