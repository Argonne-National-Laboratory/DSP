/*
 * DwMasterTrLight.cpp
 *
 *  Created on: Dec 13, 2016
 *      Author: kibaekkim
 */

#define DSP_DEBUG

#include "cplex.h"
/** Coin */
#include "OsiCpxSolverInterface.hpp"
/** Dsp */
#include <DantzigWolfe/DwMasterTrLight.h>
#include "Model/TssModel.h"

DwMasterTrLight::DwMasterTrLight(DwWorker* worker):
	DwMasterTr(worker) {}

DwMasterTrLight::~DwMasterTrLight() {}

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
#define FREE_MEMORY \
	FREE_ARRAY_PTR(delrows) \
	FREE_ARRAY_PTR(delcols) \
	FREE_ARRAY_PTR(obj_tr) \
	FREE_ARRAY_PTR(clbd_tr) \
	FREE_ARRAY_PTR(cubd_tr) \
	FREE_ARRAY_PTR(starts_tr) \
	FREE_ARRAY_PTR(rows_tr) \
	FREE_ARRAY_PTR(elements_tr) \
	FREE_ARRAY_PTR(densecol) \
	FREE_ARRAY_PTR(col_inds) \
	FREE_ARRAY_PTR(col_elems)

	int* delrows = NULL;
	int* delcols = NULL;

	/** trust-region columns */
	double* obj_tr = NULL;
	double* clbd_tr = NULL;
	double* cubd_tr = NULL;
	int* starts_tr = NULL;
	int* rows_tr = NULL;
	double* elements_tr = NULL;

	/** adding columns */
	double* densecol = NULL;
	int col_size = 0;
	int* col_inds = NULL;
	double* col_elems = NULL;

	BGN_TRY_CATCH

	if (branchobj) {
		/** remove all the branching rows */
		if (nrows_branch_ > 0) {
			delrows = new int [nrows_branch_];
			for (int i = nrows_core_; i < nrows_; ++i)
				delrows[i-nrows_core_] = i;
			si_->deleteRows(nrows_branch_, delrows);
			DSPdebugMessage("Deleted %d rows in the master.\n", nrows_branch_);
		}

		/** remove all the columns except the core TR columns */
		int ncols_tr_core = nrows_orig_*2;
		int ndelcols = si_->getNumCols() - ncols_tr_core;
		delcols = new int [ndelcols];
		for (int j = ncols_tr_core; j < si_->getNumCols(); ++j)
			delcols[j - ncols_tr_core] = j;
		si_->deleteCols(ndelcols, delcols);
		DSPdebugMessage("Deleted %d columns in the master.\n", ndelcols);

		/** count nrows_branch_ */
		nrows_branch_ = 0;
		for (unsigned j = 0; j < branchobj->index_.size(); ++j) {
			if (branchobj->lb_[j] > org_clbd_[branchobj->index_[j]]) {
				branch_row_to_col_[nrows_core_ + nrows_branch_] = branchobj->index_[j];
				si_->addRow(0, NULL, NULL, branchobj->lb_[j], COIN_DBL_MAX);
				nrows_branch_++;
			}
			if (branchobj->ub_[j] < org_cubd_[branchobj->index_[j]]) {
				branch_row_to_col_[nrows_core_ + nrows_branch_] = branchobj->index_[j];
				si_->addRow(0, NULL, NULL, -COIN_DBL_MAX, branchobj->ub_[j]);
				nrows_branch_++;
			}
		}
		DSPdebugMessage("Found %d branching rows.\n", nrows_branch_);

		/** update number of rows */
		nrows_ = nrows_core_ + nrows_branch_;

		/** change trust region */
		double* tr_center = tr_center_;
		ncols_tr_ = (nrows_orig_ + nrows_branch_) * 2;
		tr_center_ = NULL;
		tr_center_ = new double [ncols_tr_];
		CoinCopyN(tr_center, ncols_tr_core, tr_center_);
		CoinZeroN(tr_center_ + ncols_tr_core, nrows_branch_ * 2);
		FREE_ARRAY_PTR(tr_center);

		obj_tr = new double [nrows_branch_*2];
		clbd_tr = new double [nrows_branch_*2];
		cubd_tr = new double [nrows_branch_*2];
		CoinFillN(obj_tr, nrows_branch_*2, tr_size_);
		CoinZeroN(clbd_tr, nrows_branch_*2);
		CoinFillN(cubd_tr, nrows_branch_*2, COIN_DBL_MAX);

		starts_tr = new int [nrows_branch_*2 + 1];
		rows_tr = new int [nrows_branch_*2];
		elements_tr = new double [nrows_branch_*2];
		for (int j = 0, k = 0; j < nrows_branch_; ++j) {
			starts_tr[k] = k;
			rows_tr[k] = nrows_core_ + j;
			elements_tr[k] = 1.0;
			starts_tr[k+1] = k+1;
			rows_tr[k+1] = nrows_core_ + j;
			elements_tr[k+1] = -1.0;
			k += 2;
		}
		starts_tr[nrows_branch_*2] = nrows_branch_*2;
		si_->addCols(nrows_branch_*2, starts_tr, rows_tr, elements_tr, clbd_tr, cubd_tr, obj_tr);
		DSPdebugMessage("Appended %d trust-region columns in the master (%d cols).\n", nrows_branch_*2, si_->getNumCols());
		si_->writeMps("cols");

		/** add branching rows */
		for (unsigned k = 0; k < cols_generated_.size(); ++k) {
			if (cols_generated_[k]->active_) {
				col_inds = new int [nrows_];
				col_elems = new double [nrows_];
				/** get a core-row column */
				densecol = cols_generated_[k]->col_.denseVector(nrows_core_);
				/** create a column for core rows */
				col_size = 0;
				for (int i = 0; i < nrows_core_; ++i) {
					if (fabs(densecol[i]) > 1.0e-10) {
						col_inds[col_size] = i;
						col_elems[col_size] = densecol[i];
						col_size++;
					}
				}
				FREE_ARRAY_PTR(densecol);
				/** append column elements for the branching rows */
				for (unsigned j = 0, i = 0; j < branchobj->index_.size(); ++j) {
					int sparse_index = -1;
					if (branchobj->lb_[j] > org_clbd_[branchobj->index_[j]]) {
						sparse_index = cols_generated_[k]->x_.findIndex(j);
						if (sparse_index > -1) {
							double val = cols_generated_[k]->x_.getElements()[sparse_index];
							if (fabs(val) > 1.0e-10) {
								col_inds[col_size] = nrows_core_ + i;
								col_elems[col_size] = val;
								col_size++;
							}
						}
						i++;
					}
					if (branchobj->ub_[j] < org_cubd_[branchobj->index_[j]]) {
						if (sparse_index == -1)
							sparse_index = cols_generated_[k]->x_.findIndex(j);
						if (sparse_index > -1) {
							double val = cols_generated_[k]->x_.getElements()[sparse_index];
							if (fabs(val) > 1.0e-10) {
								col_inds[col_size] = nrows_core_ + i;
								col_elems[col_size] = val;
								col_size++;
							}
						}
						i++;
					}
				}
				/** assign the core-row column */
				cols_generated_[k]->col_.clear();
				cols_generated_[k]->col_.assignVector(col_size, col_inds, col_elems);
				/** add column */
				si_->addCol(cols_generated_[k]->col_, 0.0, COIN_DBL_MAX, cols_generated_[k]->obj_);
			}
		}
		DSPdebugMessage("Appended dynamic columns in the master (%d cols).\n", si_->getNumCols());
		si_->writeMps("master");

		/** restore column bounds */
		CoinCopyN(org_clbd_, ncols_orig_, node_clbd_);
		CoinCopyN(org_cubd_, ncols_orig_, node_cubd_);

		/** update column bounds at the current node */
		for (unsigned j = 0; j < branchobj->index_.size(); ++j) {
			node_clbd_[branchobj->index_[j]] = branchobj->lb_[j];
			node_cubd_[branchobj->index_[j]] = branchobj->ub_[j];
		}

		/** apply column bounds */
		worker_->setColBounds(branchobj->index_.size(),
				&branchobj->index_[0], &branchobj->lb_[0], &branchobj->ub_[0]);
	}
	END_TRY_CATCH(FREE_MEMORY)
	FREE_MEMORY
#undef FREE_MEMORY
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
	nrows_ = nrows_core_;
	nrows_branch_ = 0;

	/** allocate memory */
	obj = new double [ncols_tr_];
	clbd = new double [ncols_tr_];
	cubd = new double [ncols_tr_];
	rlbd = new double [nrows_];
	rubd = new double [nrows_];
	tr_center_ = new double [ncols_orig_];
	starts_tr = new int [ncols_tr_ + 1];
	rows_tr = new int [ncols_tr_];
	elements_tr = new double [ncols_tr_];
	node_clbd_ = new double [ncols_orig_];
	node_cubd_ = new double [ncols_orig_];
	CoinCopyN(org_clbd_, ncols_orig_, node_clbd_);
	CoinCopyN(org_cubd_, ncols_orig_, node_cubd_);

	/** initial trust region center */
	CoinZeroN(tr_center_, nrows_orig_);

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

	DSP_RTN_CHECK_RTN_CODE(initialOsiSolver());

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
