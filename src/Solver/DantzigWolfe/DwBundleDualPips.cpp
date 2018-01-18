/*
 * DwBundleDualPips.cpp
 *
 *  Created on: Dec 7, 2017
 *      Author: kibaekkim
 */

//#define DSP_DEBUG

#include "Solver/DantzigWolfe/DwBundleDualPips.h"
#include "Solver/DantzigWolfe/DwWorkerPips.h"
#include "Utility/DspUtility.h"

DwBundleDualPips& DwBundleDualPips::operator =(const DwBundleDualPips& rhs) {
	DwBundleDual::operator =(rhs);
	return *this;
}

DSP_RTN_CODE DwBundleDualPips::updateCenter(double penalty) {
	double coef;
	PipsInterface* pips = dynamic_cast<DwWorkerPips*>(worker_)->pips();

	/** update Hessian in PIPS-IPM data */
	for (int j = nrows_conv_; j < nrows_orig_; ++j) {
		/** objective coefficient */
		coef = -penalty*bestdualsol_[j];
		if (rlbd_orig_[j] > -1.0e+20)
			coef -= rlbd_orig_[j];
		if (rubd_orig_[j] < 1.0e+20)
			coef -= rubd_orig_[j];
		pips->setObjCoef(j, coef);
		//printf("PIPS ObjCoef %d %e\n", j, coef);
	}
	return DSP_RTN_OK;
}

void DwBundleDualPips::initDualSolver(
		const CoinPackedMatrix& m, 
		std::vector<double>& clbd, 
		std::vector<double>& cubd, 
		std::vector<double>& obj) {
	dynamic_cast<DwWorkerPips*>(worker_)->initPips(nrows_conv_, nrows_orig_);
}

DSP_RTN_CODE DwBundleDualPips::callMasterSolver() {

	/** run PIPS-IPM */
	dynamic_cast<DwWorkerPips*>(worker_)->solvePips(u_);

	/** TODO: parse PIPS-IPM solution status */
	status_ = DSP_STAT_OPTIMAL;

	return status_;
}

void DwBundleDualPips::assignMasterSolution(std::vector<double>& sol) {
	sol = dynamic_cast<DwWorkerPips*>(worker_)->pips()->solution_;
}

double DwBundleDualPips::getObjValue() {
	return dynamic_cast<DwWorkerPips*>(worker_)->pips()->objval_;
}

void DwBundleDualPips::addDualRow(const CoinPackedVector& v, const double lb, const double ub) {
	dynamic_cast<DwWorkerPips*>(worker_)->pips()->addRow(v, lb, ub);
}

void DwBundleDualPips::removeAllDualRows() {
	dynamic_cast<DwWorkerPips*>(worker_)->clearMats();
}

void DwBundleDualPips::printIterInfo() {
	PipsInterface* pips = dynamic_cast<DwWorkerPips*>(worker_)->pips();
	message_->print(1, "Iteration %3d: DW Bound %+e, ", itercnt_, primobj_);
	message_->print(3, "Dual %+e, Approx %+e, ", -dualobj_, -bestdualobj_-v_);
	message_->print(1, "Best Dual %+e ", -bestdualobj_);
	if (relgap_ < 1000)
		message_->print(1, "(gap %.2f %%), ", relgap_*100);
	else
		message_->print(1, "(gap Large %%), ");
	message_->print(1, "nrows %d, ncols %d, ", pips->getNumRows(), pips->getNumCols());
	message_->print(1, "timing (total %.2f, master %.2f, gencols %.2f), statue %d\n",
			t_total_, t_master_, t_colgen_, status_);
	message_->print(3, "  predicted ascent %+e, |p| %+e, alpha %+e, linerr %+e, eps %+e, u %+e, counter %d\n",
			-v_, absp_, alpha_, -linerr_, eps_, u_, counter_);

	/** log */
	log_time_.push_back(CoinGetTimeOfDay());
	log_bestdual_bounds_.push_back(-bestdualobj_);
}
