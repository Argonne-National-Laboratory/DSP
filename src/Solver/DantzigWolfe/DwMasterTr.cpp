/*
 * DwMasterTr.cpp
 *
 *  Created on: Nov 18, 2016
 *      Author: kibaekkim
 */

#define DSP_DEBUG

#include "cplex.h"
/** Coin */
#include "OsiCpxSolverInterface.hpp"
/** Dsp */
#include <DantzigWolfe/DwMasterTr.h>

DSP_RTN_CODE DwMasterTr::createProblem() {
#define FREE_MEMORY       \
	FREE_PTR(mat)

	OsiCpxSolverInterface* cpx = NULL;

	/** master problem */
	CoinPackedMatrix * mat = NULL;
	std::vector<double> obj;
	std::vector<double> clbd;
	std::vector<double> cubd;
	std::vector<double> rlbd;
	std::vector<double> rubd;

	BGN_TRY_CATCH

	/** number of initial columns */
	int ncols_tr = 2 * (nrows_orig_ + nrows_branch_);
	int ncols = cols_generated_.size();

	/** allocate memory */
	obj.reserve(ncols_tr + ncols);
	clbd.reserve(ncols_tr + ncols);
	cubd.reserve(ncols_tr + ncols);
	rlbd.reserve(nrows_branch_);
	rubd.reserve(nrows_branch_);
	rlbd_branch_ = new double [nrows_branch_];
	rubd_branch_ = new double [nrows_branch_];
	tr_center_ = new double [nrows_];

	/** initial trust region center */
	CoinZeroN(tr_center_, nrows_orig_ + nrows_branch_);
	CoinFillN(tr_center_ + nrows_orig_ + nrows_branch_, nrows_conv_, COIN_DBL_MAX);

	/** create column-wise matrix and set number of rows */
	mat = new CoinPackedMatrix(true, 0, 0);
	mat->setDimensions(nrows_, 0);

	/** add variables related to trust region */
	for (int j = 0; j < nrows_orig_ + nrows_branch_; ++j) {
		clbd.push_back(0.0);
		cubd.push_back(COIN_DBL_MAX);
		obj.push_back(tr_size_);
		double vecelem = 1.0;
		mat->appendCol(1, &j, &vecelem);

		clbd.push_back(0.0);
		cubd.push_back(COIN_DBL_MAX);
		obj.push_back(tr_size_);
		vecelem = -1.0;
		mat->appendCol(1, &j, &vecelem);
	}

	/** add initial columns */
	for (int j = 0; j < ncols; ++j) {
		clbd.push_back(cols_generated_[j]->lb_);
		cubd.push_back(cols_generated_[j]->ub_);
		obj.push_back(cols_generated_[j]->obj_);
		mat->appendCol(cols_generated_[j]->col_);
	}
	DSPdebug(mat->verifyMtx(4));

	/** Set row bounds */
	for (int i = 0; i < nrows_orig_; ++i) {
		rlbd.push_back(org_rlbd_[i]);
		rubd.push_back(org_rubd_[i]);
	}
	for (int i = 0; i < nrows_branch_; ++i) {
		int j = branch_row_to_col_[nrows_orig_ + i];
		rlbd.push_back(org_clbd_[j]);
		rubd.push_back(org_cubd_[j]);
		rlbd_branch_[i] = org_clbd_[j];
		rubd_branch_[i] = org_cubd_[j];
	}
	for (int i = 0; i < nrows_conv_; ++i) {
		rlbd.push_back(1.0);
		rubd.push_back(1.0);
	}

	/** create solver */
	si_ = new OsiCpxSolverInterface();

	/** message setting */
	//dynamic_cast<OsiClpSolverInterface*>(si_)->getModelPtr()->setLogLevel(0);
	si_->messageHandler()->logLevel(0);
	//si_->passInMessageHandler(si_->messageHandler());
	//DSPdebug(si_->messageHandler()->logLevel(4));

	/** load problem data */
	si_->loadProblem(*mat, &clbd[0], &cubd[0], &obj[0], &rlbd[0], &rubd[0]);

	/** write mps */
	DSPdebug(si_->writeMps("master"));

	/** NOTE: This does not go to Phase 1. */
	phase_ = 2;

	/** set hint parameters */
	useCpxBarrier_ = false;

	if (useCpxBarrier_) {
		cpx = dynamic_cast<OsiCpxSolverInterface*>(si_);
		/** use barrier */
		CPXsetintparam(cpx->getEnvironmentPtr(), CPX_PARAM_LPMETHOD, CPX_ALG_BARRIER);
		/** no crossover */
		CPXsetintparam(cpx->getEnvironmentPtr(), CPX_PARAM_BARCROSSALG, -1);
		CPXsetintparam(cpx->getEnvironmentPtr(), CPX_PARAM_BARDISPLAY, 0);
	} else {
		si_->setHintParam(OsiDoPresolveInResolve, true);
		si_->setHintParam(OsiDoDualInResolve, false);
		si_->setHintParam(OsiDoInBranchAndCut, false);
	}

	/** initial solve */
	//si_->initialSolve();

	/** TODO: what should be primal solution? in the original? or the restricted? */
	/** allocate memory for solution */
	bestprimsol_ = new double [ncols_orig_];
	primsol_ = new double [ncols_orig_];

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	/** release memory */
	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

DSP_RTN_CODE DwMasterTr::solve() {
#define FREE_MEMORY \
	FREE_ARRAY_PTR(piA)

	double* piA = NULL;
	int nColsAdded = 0;
	/** column generation info */
	std::vector<int> subinds;
	std::vector<int> substatus;
	std::vector<double> subcxs;
	std::vector<double> subobjs;
	std::vector<CoinPackedVector*> subsols;

	std::vector<int> delcols;

	BGN_TRY_CATCH

	/** clear master */
	if (si_->getNumCols() > 2 * (nrows_orig_ + nrows_branch_)) {
		delcols.reserve(si_->getNumCols() - 2 * (nrows_orig_ + nrows_branch_));
		for (int j = 2 * (nrows_orig_ + nrows_branch_); j < si_->getNumCols(); ++j)
			delcols.push_back(j);
		si_->deleteCols(delcols.size(), &delcols[0]);
	}
	for (unsigned i = 0; i < cols_generated_.size(); ++i)
		FREE_PTR(cols_generated_[i]);
	cols_generated_.clear();

	/** allocate memory */
	piA = new double [ncols_orig_];

	/** reset trust region */
	tr_size_ = 0.01;
	CoinZeroN(tr_center_, nrows_orig_ + nrows_branch_);
	CoinFillN(tr_center_ + nrows_orig_ + nrows_branch_, nrows_conv_, COIN_DBL_MAX);
	DSP_RTN_CHECK_RTN_CODE(updateTrustRegion());

	/** calculate pi^T A */
	DSP_RTN_CHECK_RTN_CODE(calculatePiA(tr_center_, piA));

	/** generate columns */
	DSP_RTN_CHECK_RTN_CODE(
			worker_->generateCols(phase_, piA, subinds, substatus, subcxs, subobjs, subsols));

	/** get Lagrangian dual bound */
	DSP_RTN_CHECK_RTN_CODE(getLagrangianBound(tr_center_, subobjs, dualobj_));
	DSPdebugMessage("Dual bound %e\n", dualobj_);

	/** create and add columns */
	DSP_RTN_CHECK_RTN_CODE(
			addCols(tr_center_, piA, subinds, substatus, subcxs, subobjs, subsols, nColsAdded));
	DSPdebugMessage("Number of generated columns %d\n", nColsAdded);

	/** free memory for subproblem solutions */
	for (unsigned i = 0; i < subsols.size(); ++i)
		FREE_PTR(subsols[i]);

	/** solve */
	DSP_RTN_CHECK_RTN_CODE(DwAlgo::solve());

	double relgap = (primobj_-dualobj_)/(1.0e-10+fabs(primobj_))*100;
	if (relgap > 0.01)
		status_ = DSP_STAT_MW_RESOLVE;

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

bool DwMasterTr::terminationTest(int nnewcols, int itercnt, double relgap) {

	bool term = false;
#ifdef DSP_DEBUG
	double g = 0.0;
	for (int j = 0; j < nrows_orig_ + nrows_branch_; ++j)
		g += si_->getColSolution()[j];
	DSPdebugMessage("g %e\n", g);
#endif

	if (iterlim_ <= itercnt || fabs(relgap) < 0.01) {
		status_ = DSP_STAT_FEASIBLE;
		term = true;
	}

	return term;
}

DSP_RTN_CODE DwMasterTr::updateModel(
		const double* price, /**< [in] price */
		double curLb         /**< [in] current lower bound */) {
	BGN_TRY_CATCH

	//DSPdebugMessage("current primobj %e\n", primobj_);

	if (curLb >= dualobj_ + 1.0e-4 * (primobj_ - dualobj_)) {
		message_->print(2, "  [TR] SERIOUS STEP: dual objective %e", curLb);

		if (dualobj_ > -COIN_DBL_MAX) {
			/** update proximal point */
			CoinCopyN(price, nrows_, tr_center_);
			message_->print(3, ", updated proximal point");

			/** increase trust region? */
			if (isTrBoundary(price) && curLb >= dualobj_ + 0.5 * (primobj_ - dualobj_)) {
				tr_size_ = CoinMin(2. * tr_size_, 1.0e+4);
				message_->print(3, ", increased trust region size %e", tr_size_);
			}

			/** update trust region */
			updateTrustRegion();
		}

		/** update dual bound */
		dualobj_ = curLb;
		tr_cnt_ = 0;
		message_->print(2, "\n");
	} else {
		message_->print(3, "  [TR] null step: dual objective %e", curLb);
		double rho = CoinMin(1.0, tr_size_) * (dualobj_ - curLb) / (primobj_ - dualobj_);
		if (rho > 0) tr_cnt_++;
		if (rho >= 3 || (tr_cnt_ >= 3 && fabs(rho - 2.) < 1.0))
		{
			/** decrease trust region */
			tr_size_ *= 1.0 / CoinMin(rho, 4.);
			message_->print(3, ", decreased trust region size %e", tr_size_);
			tr_cnt_ = 0;

			/** update trust region */
			updateTrustRegion();
		}
		message_->print(3, "\n");
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

bool DwMasterTr::isTrBoundary(const double* price) {
	bool isBoundary = false;

	BGN_TRY_CATCH

	for (int j = 0; j < nrows_orig_ + nrows_branch_; ++j)
		if (fabs(price[j] - tr_center_[j]) < 1.0e-10) {
			isBoundary = true;
			break;
		}

	END_TRY_CATCH_RTN(;,isBoundary)

	return isBoundary;
}

DSP_RTN_CODE DwMasterTr::updateTrustRegion() {
	BGN_TRY_CATCH

	for (int j = 0; j < nrows_orig_ + nrows_branch_; ++j) {
		si_->setObjCoeff(j*2, tr_center_[j] + tr_size_);
		si_->setObjCoeff(j*2+1, -tr_center_[j] + tr_size_);
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}
