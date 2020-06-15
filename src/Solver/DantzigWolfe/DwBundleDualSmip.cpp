/*
 * DwBundleDualSmip.cpp
 *
 *  Created on: Jun 27, 2018
 *      Author: Kibaek Kim
 */

//#define DSP_DEBUG
#include "SolverInterface/DspOsiClp.h"
#include "SolverInterface/DspOsiCpx.h"
#include "SolverInterface/DspOsiGrb.h"
#include "Solver/DantzigWolfe/DwBundleDualSmip.h"
#include "Utility/DspUtility.h"

DwBundleDualSmip::DwBundleDualSmip(DwWorker* worker):
	DwBundleDual(worker) {}

DwBundleDualSmip::DwBundleDualSmip(const DwBundleDualSmip& rhs):
DwBundleDual(rhs), tss_(rhs.tss_) {}

DwBundleDualSmip& DwBundleDualSmip::operator =(const DwBundleDualSmip& rhs) {
	DwBundleDual::operator =(rhs);
	tss_ = rhs.tss_;
	return *this;
}

DwBundleDualSmip::~DwBundleDualSmip() {
	tss_ = NULL;
}

DSP_RTN_CODE DwBundleDualSmip::init() {

	BGN_TRY_CATCH

	/** retrieve tss model */
	tss_ = dynamic_cast<TssModel*>(model_);

	DSP_RTN_CHECK_THROW(DwMaster::init());

	/** modify non-anticipativity constraints in the following form:
	   -x + x_1        = 0
	   -x +    x_2     = 0
	   -x +        x_3 = 0
	*/
	std::vector<int> delrows(mat_orig_->getNumRows(), 0);
	std::iota(delrows.begin(), delrows.end(), 0);
	mat_orig_->deleteRows(delrows.size(), &delrows[0]);

	int indices[1];
	double elements[] = {1.0};
	for (int i = 0; i < tss_->getNumScenarios(); ++i) {
		for (int j = 0; j < tss_->getNumCols(0); ++j) {
			indices[0] = i * tss_->getNumCols(0) + j;
			mat_orig_->appendRow(1, indices, elements);
		}
	}
	DSPdebug(mat_orig_->verifyMtx(4));

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DwBundleDualSmip::createPrimalProblem() {
	BGN_TRY_CATCH

	/** create column-wise matrix and set number of rows */
	std::shared_ptr<CoinPackedMatrix> mat(new CoinPackedMatrix(true, 0, 0));
	mat->setDimensions(nrows_, 0);
	std::vector<int> indices(tss_->getNumScenarios(), 0);
	std::vector<double> elements(tss_->getNumScenarios(), -1.0);
	for (int j = 0; j < tss_->getNumCols(0); ++j) {
		for (int s = 0; s < tss_->getNumScenarios(); ++s)
			indices[s] = tss_->getNumScenarios() + s * tss_->getNumCols(0) + j;
		mat->appendCol(tss_->getNumScenarios(), &indices[0], &elements[0]);
	}

	/** row bounds */
	std::vector<double> clbd(tss_->getNumCols(0), -COIN_DBL_MAX);
	std::vector<double> cubd(tss_->getNumCols(0), +COIN_DBL_MAX);
	std::vector<double> obj(tss_->getNumCols(0), 0.0);
	std::vector<double> rlbd(nrows_, 0.0);
	std::vector<double> rubd(nrows_, 0.0);
	std::fill(rlbd.begin(), rlbd.begin() + nrows_conv_, 1.0);
	std::fill(rubd.begin(), rubd.begin() + nrows_conv_, 1.0);

	/** create solver */
	DspOsi * osi = createDspOsi();
	if (osi)
		primal_si_.reset(osi);
	else
		throw CoinError("Failed to create DspOsi", "createPrimalProblem", "DwBundleDualSmip");

	/** load problem data */
	primal_si_->si_->loadProblem(*mat, &clbd[0], &cubd[0], &obj[0], &rlbd[0], &rubd[0]);

	/** set display */
	primal_si_->setLogLevel(par_->getIntParam("DW/MASTER/SOLVER/LOG_LEVEL"));

	/** set number of cores */
	primal_si_->setNumCores(par_->getIntParam("NUM_CORES"));

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DwBundleDualSmip::createDualProblem() {
	BGN_TRY_CATCH

	/** master problem */
	std::shared_ptr<CoinPackedMatrix> mat(nullptr);
	std::vector<double> clbd(nrows_), cubd(nrows_), obj(nrows_);

	/** necessary dual variables */
	dualsol_.resize(nrows_);
	std::fill(dualsol_.begin(), dualsol_.begin() + nrows_conv_, COIN_DBL_MAX);
	std::fill(dualsol_.begin() + nrows_conv_, dualsol_.end(), 0.0);
	bestdualsol_ = dualsol_;

	/** other initialization */
	p_.reserve(nrows_orig_+nrows_branch_);
	d_.reserve(nrows_orig_+nrows_branch_);
	std::fill(p_.begin(), p_.end(), 0.0);
	std::fill(d_.begin(), d_.end(), 0.0);

	/** create row-wise matrix and set number of rows */
	mat.reset(new CoinPackedMatrix(false, 0, 0));
	mat->setDimensions(0, nrows_);

	std::vector<int> indices;
	std::vector<double> elements(tss_->getNumScenarios(), 1.0);
	indices.reserve(tss_->getNumScenarios());
	for (int j = 0; j < tss_->getNumCols(0); ++j) {
		indices.clear();
		for (int s = 0; s < tss_->getNumScenarios(); ++s)
			indices.push_back(tss_->getNumScenarios() + s * tss_->getNumCols(0) + j);
		mat->appendRow(indices.size(), &indices[0], &elements[0]);
	}
	numFixedRows_ = tss_->getNumCols(0);

	std::vector<double> rlbd(tss_->getNumCols(0), 0.0);
	std::vector<double> rubd(tss_->getNumCols(0), 0.0);

	std::fill(clbd.begin(), clbd.end(), -COIN_DBL_MAX);
	std::fill(cubd.begin(), cubd.end(), +COIN_DBL_MAX);
	std::fill(obj.begin(), obj.begin() + nrows_conv_, -1.0);
	for (int i = 0; i < nrows_orig_; ++i) {
		obj[nrows_conv_+i] = -u_*bestdualsol_[nrows_conv_+i];
	}

	/** initialize external solver */
	initDualSolver(*mat, clbd, cubd, obj, rlbd, rubd);

	/** set quadratic objective term */
	DSP_RTN_CHECK_THROW(updateCenter(u_));

	/** set number of cores */
	primal_si_->setNumCores(par_->getIntParam("NUM_CORES"));

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

void DwBundleDualSmip::removeAllPrimCols() {
	std::vector<int> delcols(primal_si_->si_->getNumCols() - tss_->getNumCols(0));
	std::iota(delcols.begin(), delcols.end(), tss_->getNumCols(0));
	primal_si_->si_->deleteCols(primal_si_->si_->getNumCols() - tss_->getNumCols(0), &delcols[0]);
}

void DwBundleDualSmip::removeAllDualRows() {
	std::vector<int> delrows(getSiPtr()->getNumRows() - tss_->getNumCols(0));
	std::iota(delrows.begin(), delrows.end(), tss_->getNumCols(0));
	getSiPtr()->deleteRows(getSiPtr()->getNumRows() - tss_->getNumCols(0), &delrows[0]);
}
