/*
 * DwMasterTrLight.cpp
 *
 *  Created on: Dec 13, 2016
 *      Author: kibaekkim
 */

//#define DSP_DEBUG

#include "cplex.h"
/** Coin */
#include "OsiCpxSolverInterface.hpp"
/** Dsp */
#include <DantzigWolfe/DwMasterTrLight.h>
#include "Model/TssModel.h"

DwMasterTrLight::DwMasterTrLight(DwWorker* worker):
DwMasterTr(worker),
node_clbd_(NULL),
node_cubd_(NULL) {}

DwMasterTrLight::~DwMasterTrLight() {
	FREE_ARRAY_PTR(node_clbd_);
	FREE_ARRAY_PTR(node_cubd_);
}

bool DwMasterTrLight::chooseBranchingObjects(
		DspBranch*& branchingUp, /**< [out] branching-up object */
		DspBranch*& branchingDn  /**< [out] branching-down object */) {
	int findPhase = 0;
	bool branched = false;
	double dist, maxdist = 1.0e-6;
	int branchingIndex = -1;
	double branchingValue;

	BGN_TRY_CATCH

	/** cleanup */
	FREE_PTR(branchingUp)
	FREE_PTR(branchingDn)

	if (model_->isStochastic()) {
		TssModel* tss = dynamic_cast<TssModel*>(model_);

		/** most fractional value */
		for (int j = 0; j < tss->getNumScenarios() * tss->getNumCols(0); ++j) {
			if (org_ctype_[j] == 'C') continue;
			dist = fabs(primsol_[j] - floor(primsol_[j] + 0.5));
			if (dist > maxdist) {
				maxdist = dist;
				branchingIndex = j;
				branchingValue = primsol_[j];
			}
		}

		if (branchingIndex < 0) {
			/** most fractional value */
			for (int j = tss->getNumScenarios() * tss->getNumCols(0); j < ncols_orig_; ++j) {
				if (org_ctype_[j] == 'C') continue;
				dist = fabs(primsol_[j] - floor(primsol_[j] + 0.5));
				if (dist > maxdist) {
					maxdist = dist;
					branchingIndex = j;
					branchingValue = primsol_[j];
				}
			}
		}
	} else {
		findPhase = 0;
		while (findPhase < 2 && branchingIndex < 0) {
			/** most fractional value */
			for (int j = 0; j < ncols_orig_; ++j) {
				if (org_ctype_[j] == 'C') continue;
				/** do not consider those with non-negative objective coefficients */
				if (findPhase == 0 && org_obj_[j] >= 0) continue;
				/** do not consider those with negative objective coefficients */
				if (findPhase == 1 && org_obj_[j] < 0) continue;
				dist = fabs(primsol_[j] - floor(primsol_[j] + 0.5));
				if (dist > maxdist) {
					maxdist = dist;
					branchingIndex = j;
					branchingValue = primsol_[j];
				}
			}
			findPhase++;
		}
	}

	if (branchingIndex > -1) {
		DSPdebugMessage("Creating branch objects on column %d (value %e).\n", branchingIndex, branchingValue);
		branched = true;

		/** creating branching objects */
		branchingUp = new DspBranch();
		branchingDn = new DspBranch();
		for (int j = 0; j < ncols_orig_; ++j) {
			if (org_ctype_[j] == 'C') continue;
			if (branchingIndex == j) {
				branchingUp->push_back(j, ceil(branchingValue), node_cubd_[j]);
				branchingDn->push_back(j, node_clbd_[j], floor(branchingValue));
			} else if (node_clbd_[j] > org_clbd_[j] || node_cubd_[j] < org_cubd_[j]) {
				/** store any bound changes made in parent nodes */
				branchingUp->push_back(j, node_clbd_[j], node_cubd_[j]);
				branchingDn->push_back(j, node_clbd_[j], node_cubd_[j]);
			}
		}
	} else {
		DSPdebugMessage("No branch object is found.\n");
		//DSPdebug(si_->writeMps("Incumbent"));
	}

	END_TRY_CATCH_RTN(;,false)

	return branched;
}

void DwMasterTrLight::setBranchingObjects(const DspBranch* branchobj) {
	BGN_TRY_CATCH
	if (branchobj) {
		/** restore column bounds */
		for (int j = 0; j < ncols_orig_; ++j) {
			if (org_ctype_[j] == 'C') continue;
			node_clbd_[j] = org_clbd_[j];
			node_cubd_[j] = org_cubd_[j];
		}

		/** update column bounds at the current node */
		for (unsigned j = 0; j < branchobj->index_.size(); ++j) {
			node_clbd_[branchobj->index_[j]] = branchobj->lb_[j];
			node_cubd_[branchobj->index_[j]] = branchobj->ub_[j];
		}

		/** apply column bounds */
		worker_->setColBounds(branchobj->index_.size(),
				&branchobj->index_[0], &branchobj->lb_[0], &branchobj->ub_[0]);

		/** filter generated columns */
		for (unsigned k = 0; k < cols_generated_.size(); ++k) {
			for (int j = 0; j < ncols_orig_; ++j) {
				if (org_ctype_[j] == 'C') continue;
				double xval = cols_generated_[k]->x_.isExistingIndex(j) ? cols_generated_[k]->x_[j] : 0.0;
				if (xval < node_clbd_[j] || xval > node_cubd_[j])
					cols_generated_[k]->active_ = false;
				else
					cols_generated_[k]->active_ = true;
			}
		}
	}
	END_TRY_CATCH(;)
}

DSP_RTN_CODE DwMasterTrLight::createProblem() {
#define FREE_MEMORY \
	FREE_PTR(mat) \
	FREE_ARRAY_PTR(obj) \
	FREE_ARRAY_PTR(clbd) \
	FREE_ARRAY_PTR(cubd) \
	FREE_ARRAY_PTR(rlbd) \
	FREE_ARRAY_PTR(rubd) \
	FREE_ARRAY_PTR(starts_tr) \
	FREE_ARRAY_PTR(rows_tr) \
	FREE_ARRAY_PTR(elements_tr)

	OsiCpxSolverInterface* cpx = NULL;

	/** master problem */
	CoinPackedMatrix * mat = NULL;
	double* obj = NULL;
	double* clbd = NULL;
	double* cubd = NULL;
	double* rlbd = NULL;
	double* rubd = NULL;

	/** trust-region columns */
	int* starts_tr = NULL;
	int* rows_tr = NULL;
	double* elements_tr = NULL;

	BGN_TRY_CATCH

	/** number of initial columns */
	ncols_tr_ = 2 * nrows_orig_;

	/** column index of the first generated columns added */
	ncols_start_ = ncols_tr_;

	/** do add consider branch rows */
	nrows_ -= nrows_branch_;
	nrows_branch_ = 0;

	/** allocate memory */
	obj = new double [ncols_tr_];
	clbd = new double [ncols_tr_];
	cubd = new double [ncols_tr_];
	rlbd = new double [nrows_];
	rubd = new double [nrows_];
	tr_center_ = new double [nrows_];
	starts_tr = new int [ncols_tr_ + 1];
	rows_tr = new int [ncols_tr_];
	elements_tr = new double [ncols_tr_];
	node_clbd_ = new double [ncols_orig_];
	node_cubd_ = new double [ncols_orig_];
	CoinCopyN(org_clbd_, ncols_orig_, node_clbd_);
	CoinCopyN(org_cubd_, ncols_orig_, node_cubd_);

	/** initial trust region center */
	CoinZeroN(tr_center_, nrows_orig_);
	CoinFillN(tr_center_ + nrows_orig_, nrows_conv_, COIN_DBL_MAX);

	/** create column-wise matrix and set number of rows */
	mat = new CoinPackedMatrix(true, 0, 0);
	mat->setDimensions(nrows_, 0);

	/** add variables related to trust region */
	CoinFillN(obj, ncols_tr_, tr_size_);
	CoinZeroN(clbd, ncols_tr_);
	CoinFillN(cubd, ncols_tr_, COIN_DBL_MAX);
	for (int j = 0, k = 0; j < nrows_orig_; ++j) {
		starts_tr[k] = k;
		rows_tr[k] = j;
		elements_tr[k] = 1.0;
		starts_tr[k+1] = k+1;
		rows_tr[k+1] = j;
		elements_tr[k+1] = -1.0;
		k += 2;
	}
	starts_tr[ncols_tr_] = ncols_tr_;
	mat->appendCols(ncols_tr_, starts_tr, rows_tr, elements_tr);

	DSPdebug(mat->verifyMtx(4));

	/** Set row bounds */
	CoinCopyN(org_rlbd_, nrows_orig_, rlbd);
	CoinCopyN(org_rubd_, nrows_orig_, rubd);
	CoinFillN(rlbd + nrows_orig_, nrows_conv_, 1.0);
	CoinFillN(rubd + nrows_orig_, nrows_conv_, 1.0);

	/** create solver */
	si_ = new OsiCpxSolverInterface();

	/** message setting */
	si_->messageHandler()->logLevel(0);
	DSPdebug(si_->messageHandler()->logLevel(4));

	/** load problem data */
	si_->loadProblem(*mat, clbd, cubd, obj, rlbd, rubd);

	/** write mps */
	DSPdebug(si_->writeMps("master"));

	/** NOTE: This does not go to Phase 1. */
	phase_ = 2;

	/** set hint parameters */
	useCpxBarrier_ = par_->getBoolParam("DW/MASTER/IPM");

	cpx = dynamic_cast<OsiCpxSolverInterface*>(si_);
	if (cpx)
		CPXsetintparam(cpx->getEnvironmentPtr(), CPX_PARAM_THREADS, 1);

	if (useCpxBarrier_) {
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
