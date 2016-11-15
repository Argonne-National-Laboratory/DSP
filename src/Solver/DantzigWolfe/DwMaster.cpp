//
// Created by Kibaek Kim on 8/27/16.
//

#define DSP_DEBUG

#include "cplex.h"
/** Coin */
#include "OsiClpSolverInterface.hpp"
#include "OsiCpxSolverInterface.hpp"
/** Dsp */
#include "Utility/DspUtility.h"
#include "Solver/DantzigWolfe/DwMaster.h"
#include "Model/TssModel.h"

DSP_RTN_CODE DwMaster::init() {
#define FREE_MEMORY \
	FREE_PTR(org_mat); \
	FREE_ARRAY_PTR(org_clbd); \
	FREE_ARRAY_PTR(org_cubd); \
	FREE_ARRAY_PTR(org_obj); \
	FREE_ARRAY_PTR(org_ctype); \
	FREE_ARRAY_PTR(org_rlbd); \
	FREE_ARRAY_PTR(org_rubd);

	TssModel* tss = NULL;
	CoinPackedMatrix* org_mat = NULL;
	double* org_clbd = NULL;
	double* org_cubd = NULL;
	double* org_obj = NULL;
	char* org_ctype = NULL;
	double* org_rlbd = NULL;
	double* org_rubd = NULL;

	BGN_TRY_CATCH

	/** feasibility pump solver */
	fpump_ = new DwFeasPump(model_, par_, message_);

	/** subproblem solver */
    worker_ = new DwWorker(model_, par_, message_);

    if (model_->isStochastic()) {
    	DSPdebugMessage("Loading stochastic model.\n");

    	/** two-stage stochastic model */
    	tss = dynamic_cast<TssModel*>(model_);

    	/** get DE model */
    	DSP_RTN_CHECK_THROW(model_->getFullModel(org_mat, org_clbd, org_cubd, org_ctype, org_obj, org_rlbd, org_rubd));

    	int nscen = tss->getNumScenarios();
    	int ncols_first_stage = tss->getNumCols(0);
    	int ncols = org_mat->getNumCols() + ncols_first_stage * (nscen - 1);
    	const double* probability = tss->getProbability();

    	org_mat_ = new CoinPackedMatrix(org_mat->isColOrdered(), 0, 0);
    	org_mat_->setDimensions(0, ncols);

    	/** add non-anticipativity constraints */
    	int indices[2];
    	double elements[] = {1.0, -1.0};
    	for (int i = 0; i < nscen; ++i) {
    		if (i < nscen - 1) {
        		for (int j = 0; j < ncols_first_stage; ++j) {
        			indices[0] = i * ncols_first_stage + j;
        			indices[1] = (i+1) * ncols_first_stage + j;
            		org_mat_->appendRow(2, indices, elements);
        		}
    		} else {
        		for (int j = 0; j < ncols_first_stage; ++j) {
        			indices[0] = i * ncols_first_stage + j;
        			indices[1] = j;
            		org_mat_->appendRow(2, indices, elements);
        		}
    		}
    	}
    	DSPdebug(org_mat_->verifyMtx(4));

    	org_clbd_ = new double [ncols];
    	org_cubd_ = new double [ncols];
    	org_ctype_ = new char [ncols];
    	org_obj_ = new double [ncols];
    	org_rlbd_ = new double [org_mat_->getNumRows()];
    	org_rubd_ = new double [org_mat_->getNumRows()];
		for (int s = 0; s < nscen; ++s) {
			CoinCopyN(org_clbd, ncols_first_stage, org_clbd_ + s * ncols_first_stage);
			CoinCopyN(org_cubd, ncols_first_stage, org_cubd_ + s * ncols_first_stage);
			CoinCopyN(org_ctype, ncols_first_stage, org_ctype_ + s * ncols_first_stage);
	    	for (int j = 0; j < ncols_first_stage; ++j)
	    		org_obj_[s * ncols_first_stage + j] = org_obj[j] * probability[s];
		}
		CoinCopyN(org_clbd + ncols_first_stage,  ncols - nscen * ncols_first_stage,
				org_clbd_ + nscen * ncols_first_stage);
		CoinCopyN(org_cubd + ncols_first_stage,  ncols - nscen * ncols_first_stage,
				org_cubd_ + nscen * ncols_first_stage);
		CoinCopyN(org_ctype + ncols_first_stage, ncols - nscen * ncols_first_stage,
				org_ctype_ + nscen * ncols_first_stage);
		CoinZeroN(org_obj_ + nscen * ncols_first_stage, ncols - nscen * ncols_first_stage);
		CoinZeroN(org_rlbd_, org_mat_->getNumRows());
		CoinZeroN(org_rubd_, org_mat_->getNumRows());
    } else {
    	/** retrieve the original master problem structure */
    	model_->decompose(0, NULL, 0, NULL, NULL, NULL,
    			org_mat_, org_clbd_, org_cubd_, org_ctype_, org_obj_, org_rlbd_, org_rubd_);
    }

	ncols_orig_ = org_mat_->getNumCols(); /**< number of columns in the original master */
	nrows_orig_ = org_mat_->getNumRows(); /**< number of rows in the original master */
	nrows_branch_ = 0; /**< number of branching rows in the restricted master */
	nrows_conv_ = model_->getNumSubproblems(); /**< number of convex combination rows in the restricted master */

	/** maps each branching row to original column index */
	for (int j = 0; j < ncols_orig_; ++j)
		if (org_ctype_[j] != 'C') {
			branch_row_to_col_[nrows_orig_ + nrows_branch_] = j;
			nrows_branch_++;
		}

	/** number of rows in the restricted master */
	nrows_ = nrows_orig_ + nrows_branch_ + nrows_conv_;

	DSPdebugMessage("nrwos_ %d, nrows_orig_ %d, nrows_branch_ %d, nrows_conv_ %d\n",
			nrows_, nrows_orig_, nrows_branch_, nrows_conv_);

#if 0
	/** maps each subproblem to branching constraints */
	if (model_->isStochastic()) {
		for (int s = 0; s < tss->getNumScenarios(); ++s) {
			std::vector<int> ccols;
			/** first-stage variables */
			for (int j = 0; j < tss->getNumCols(0); ++j) {
				int ccol = s * tss->getNumCols(0) + j;
				if (org_ctype_[ccol] != 'C')
					ccols.push_back(ccol);
			}
			/** second-stage variables */
			for (int j = 0; j < tss->getNumCols(1); ++j) {
				int ccol = tss->getNumScenarios() * tss->getNumCols(0) + s * tss->getNumCols(1) + j;
				if (org_ctype_[ccol] != 'C')
					ccols.push_back(ccol);
			}
			subproblem_to_branch_rows_[s] = ccols;
		}
	} else {
		for (int s = 0; s < model_->getNumSubproblems(); ++s) {
			std::vector<int> ccols;
			for (int j = 0; j < model_->getNumSubproblemCouplingCols(s); ++j) {
				int ccol = model_->getSubproblemCouplingColIndices(s)[j];
				if (org_ctype_[ccol] != 'C')
					ccols.push_back(ccol);
			}
			subproblem_to_branch_rows_[s] = ccols;
		}
	}
#endif

	/** generate initial columns */
	DSP_RTN_CHECK_RTN_CODE(initialColumns());

	/** create problem */
	DSP_RTN_CHECK_RTN_CODE(createProblem());

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

DSP_RTN_CODE DwMaster::createProblem() {
#define FREE_MEMORY       \
	FREE_PTR(mat)         \
	FREE_ARRAY_PTR(obj)   \
	FREE_ARRAY_PTR(clbd)  \
	FREE_ARRAY_PTR(cubd)  \
	FREE_ARRAY_PTR(rlbd)  \
	FREE_ARRAY_PTR(rubd)

	/** master problem */
	CoinPackedMatrix * mat = NULL;
	double * obj = NULL;
	double * clbd = NULL;
	double * cubd = NULL;
	double * rlbd = NULL;
	double * rubd = NULL;

	BGN_TRY_CATCH

	/** number of initial columns */
	int ncols = cols_generated_.size();

	/** allocate memory */
	clbd = new double [ncols];
	cubd = new double [ncols];
	obj  = new double [ncols];
	rlbd = new double [nrows_];
	rubd = new double [nrows_];
	rlbd_branch_ = new double [nrows_branch_];
	rubd_branch_ = new double [nrows_branch_];

	/** create column-wise matrix and set number of rows */
	mat = new CoinPackedMatrix(true, 0, 0);
	mat->setDimensions(nrows_, 0);

	/** add initial columns */
	for (int j = 0; j < ncols; ++j) {
		clbd[j] = cols_generated_[j]->lb_;
		cubd[j] = cols_generated_[j]->ub_;
		obj[j] = cols_generated_[j]->obj_;
		mat->appendCol(cols_generated_[j]->col_);
	}
	DSPdebug(mat->verifyMtx(4));

	/** Set row bounds */
	CoinCopyN(org_rlbd_, nrows_orig_, rlbd);
	CoinCopyN(org_rubd_, nrows_orig_, rubd);
	for (int i = 0; i < nrows_branch_; ++i) {
		int j = branch_row_to_col_[nrows_orig_ + i];
		rlbd[nrows_orig_+i] = org_clbd_[j];
		rubd[nrows_orig_+i] = org_cubd_[j];
		rlbd_branch_[i] = org_clbd_[j];
		rubd_branch_[i] = org_cubd_[j];
	}
	CoinFillN(rlbd + nrows_orig_ + nrows_branch_, nrows_conv_, 1.0);
	CoinFillN(rubd + nrows_orig_ + nrows_branch_, nrows_conv_, 1.0);

	/** create solver */
	si_ = new OsiCpxSolverInterface();
	//dynamic_cast<OsiClpSolverInterface*>(si_)->getModelPtr()->setLogLevel(0);
	si_->messageHandler()->logLevel(0);
	//DSPdebug(si_->messageHandler()->logLevel(4));

	/** load problem data */
	si_->loadProblem(*mat, clbd, cubd, obj, rlbd, rubd);

	/** create auxiliary columns for Phase 1 problem */
	int auxcolind;
	double auxcolval;
	for (int i = 0; i < nrows_orig_ + nrows_branch_; ++i) {
		auxcolind = i;
		switch (si_->getRowSense()[i]) {
		case 'G':
			auxcolval = 1.0;
			auxcols_.push_back(new CoinPackedVector(1, &auxcolind, &auxcolval));
			break;
		case 'L':
			auxcolval = -1.0;
			auxcols_.push_back(new CoinPackedVector(1, &auxcolind, &auxcolval));
			break;
		case 'E':
		case 'R':
			auxcolval = 1.0;
			auxcols_.push_back(new CoinPackedVector(1, &auxcolind, &auxcolval));
			auxcolval = -1.0;
			auxcols_.push_back(new CoinPackedVector(1, &auxcolind, &auxcolval));
			break;
		default:
			break;
		}
	}
	phase_ = 2;

	/** write mps */
	DSPdebug(si_->writeMps("master"));

	/** set hint parameters */
	si_->setHintParam(OsiDoPresolveInResolve, true);
	si_->setHintParam(OsiDoDualInResolve, false);
	si_->setHintParam(OsiDoInBranchAndCut, false);

	/** initial solve */
	si_->initialSolve();

	if (si_->isProvenPrimalInfeasible()) {
		DSPdebugMessage("The initial master is infeasible. Phase I problem is solved.\n");
		/** set phase 1 problem */
		for (int j = 0; j < si_->getNumCols(); ++j)
			si_->setObjCoeff(j, 0.0);
		auxcolindices_.clear();
		for (unsigned j = 0; j < auxcols_.size(); ++j) {
			si_->addCol(*auxcols_[j], 0.0, COIN_DBL_MAX, 1.0);
			auxcolindices_.push_back(ncols + j);
		}
		phase_ = 1;
		DSPdebug(si_->writeMps("initialPhase1"));

		/** reoptimize */
		si_->resolve();
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

DSP_RTN_CODE DwMaster::heuristics() {
	return DSP_RTN_OK;
	BGN_TRY_CATCH

	if (phase_ == 2 && runHeuristics_) {
		/** run feasibility pump */
		DSP_RTN_CHECK_RTN_CODE(fpump_->copyColumns(cols_generated_));
		DSP_RTN_CHECK_RTN_CODE(fpump_->setBranchRowBounds(
				si_->getRowLower() + nrows_orig_, si_->getRowUpper() + nrows_orig_));
		fpump_->setBestPrimalObjective(bestprimobj_);
		DSP_RTN_CHECK_RTN_CODE(fpump_->solve());
		if (fpump_->getBestPrimalObjective() < bestprimobj_ - 1.0e-10) {
			bestprimobj_ = fpump_->getBestPrimalObjective();
			CoinCopyN(fpump_->getBestPrimalSolution(), ncols_orig_, bestprimsol_);
		}
		int ncols_before = cols_generated_.size();
		DSP_RTN_CHECK_RTN_CODE(fpump_->getNewCols(cols_generated_));
		if (ncols_before < cols_generated_.size()) {
			for (int i = ncols_before; i < cols_generated_.size(); ++i)
				if (cols_generated_[i]->active_) {
					si_->addCol(cols_generated_[i]->col_,
							cols_generated_[i]->lb_, cols_generated_[i]->ub_, cols_generated_[i]->obj_);
				}
			status_ = DSP_STAT_MW_RESOLVE;
		}
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DwMaster::hDiving() {
#define FREE_MEMORY \
	FREE_ARRAY_PTR(rlbd) \
	FREE_ARRAY_PTR(rubd)

	/** original row bounds for branching constraints */
	double* rlbd = NULL;
	double* rubd = NULL;

	double inttol = 1.0e-8; /**< integer tolerance */
	double minfrac = 0.5; /**< minimum fractional value */
	int minfracind; /**< minimum fractional indices */
	double newlbd, newubd; /**< new bounds */

	BGN_TRY_CATCH

	/** copy bounds */
	rlbd = new double [nrows_branch_];
	rubd = new double [nrows_branch_];
	CoinCopyN(si_->getRowLower() + nrows_orig_, nrows_branch_, rlbd);
	CoinCopyN(si_->getRowUpper() + nrows_orig_, nrows_branch_, rubd);

	/** find lowest fractional variable */
	minfracind = -1;
	for (int i = 0; i < nrows_branch_; ++i) {
		double x = si_->getRowActivity()[nrows_orig_+i];
		double frac = fabs(x - floor(x + 0.5));
		if (frac < minfrac && frac > inttol) {
			minfrac = frac;
			minfracind = i;
			if (frac - floor(frac) < 0.5) {
				newlbd = si_->getRowLower()[nrows_orig_+i];
				newubd = floor(frac);
			} else {
				newlbd = ceil(frac);
				newubd = si_->getRowUpper()[nrows_orig_+i];
			}
		}
	}

	if (minfracind > -1) {
		/** bound variable */
		si_->setRowBounds(nrows_orig_ + minfracind, newlbd, newubd);
		worker_->setColBounds(branch_row_to_col_[nrows_orig_+minfracind], newlbd, newubd);

		/** resolve */
	} else {
		/** found integer variable */
	}

	/** restore bounds */
	for (int i = 0; i < nrows_branch_; ++i) {
		si_->setRowBounds(nrows_orig_ + i, rlbd[i], rubd[i]);
		worker_->setColBounds(branch_row_to_col_[nrows_orig_+i], rlbd[i], rubd[i]);
	}

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

DSP_RTN_CODE DwMaster::solveMip() {
#define FREE_MEMORY \
	FREE_PTR(si)

	OsiSolverInterface* si = NULL;

	BGN_TRY_CATCH

	/** copy si */
	si = si_->clone(true);

	/** set solution */
	si->setColSolution(si_->getColSolution());

	int ncols = si->getNumCols();
	for (int i = 0; i < nrows_branch_; ++i) {
		int row = nrows_orig_ + i;
		double elem = -1.0;
		/** add auxiliary integer variables */
		si->addCol(1, &row, &elem,
				si->getRowLower()[row],
				si->getRowUpper()[row], 0.0);
		/** modify row bounds for branching constraints */
		si->setRowBounds(row, 0.0, 0.0);
		/** set integer */
		si->setInteger(ncols+i);
	}

	OsiCpxSolverInterface*cpx = dynamic_cast<OsiCpxSolverInterface*>(si);
	if (cpx) cpx->switchToMIP();

	DSPdebug(si->writeMps("MasterMip"));

	/** solve */
	si->branchAndBound();

	if (si->isProvenOptimal()) {
		bestprimobj_ = si->getObjValue();
		/** recover original solution */
		CoinZeroN(bestprimsol_, ncols_orig_);
		for (unsigned k = 0, j = 0; k < cols_generated_.size(); ++k) {
			/** do not consider inactive columns */
			if (cols_generated_[k]->active_ == false)
				continue;
			CoinPackedVector xlam = cols_generated_[k]->x_ * si->getColSolution()[j];
			for (int i = 0; i < xlam.getNumElements(); ++i)
				bestprimsol_[xlam.getIndices()[i]] += xlam.getElements()[i];
			j++;
		}
	}

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

bool DwMaster::chooseBranchingObjects(
		DspBranch*& branchingUp, /**< [out] branching-up object */
		DspBranch*& branchingDn  /**< [out] branching-down object */) {

	bool branched = false;
	double x, dist, maxdist = 1.0e-6;
	int branchingIndex = -1;
	double branchingValue;

	BGN_TRY_CATCH

	const double* rlbd = si_->getRowLower();
	const double* rubd = si_->getRowUpper();
	const double* Ax = si_->getRowActivity();

	/** cleanup */
	FREE_PTR(branchingUp)
	FREE_PTR(branchingDn)

	/** most fractional value */
	for (unsigned j = 0; j < nrows_branch_; ++j) {
		x = Ax[nrows_orig_ + j];
		//DSPdebugMessage("Checking integrality: row %d, val %e\n", nrows_orig_ + j, x);
		dist = fabs(x - floor(x + 0.5));
		if (dist > maxdist) {
			maxdist = dist;
			branchingIndex = nrows_orig_ + j;
			branchingValue = x;
		}
	}

	if (branchingIndex > -1) {
		DSPdebugMessage("Creating branch objects on column %d (value %e).\n", branchingIndex - nrows_orig_, branchingValue);
		branched = true;

		/** creating branching objects */
		branchingUp = new DspBranch();
		branchingDn = new DspBranch();
		for (unsigned j = nrows_orig_; j < nrows_orig_ + nrows_branch_; ++j) {
			if (branchingIndex == j) {
				branchingUp->push_back(branchingIndex, ceil(branchingValue), rubd[branchingIndex]);
				branchingDn->push_back(branchingIndex, rlbd[branchingIndex], floor(branchingValue));
			} else if (rlbd[j] > rlbd_branch_[j - nrows_orig_] || rubd[j] < rubd_branch_[j - nrows_orig_]) {
				/** store any bound changes made in parent nodes */
				branchingUp->push_back(j, rlbd[j], rubd[j]);
				branchingDn->push_back(j, rlbd[j], rubd[j]);
			}
		}
	} else {
		//DSPdebug(si_->writeMps("Incumbent"));
	}

	END_TRY_CATCH_RTN(;,false)

	return branched;
}

void DwMaster::setBranchingObjects(const DspBranch* branchobj) {
	if (branchobj) {
		/** restore original bounds */
		for (int j = 0; j < nrows_branch_; ++j) {
			si_->setRowBounds(nrows_orig_ + j, rlbd_branch_[j], rubd_branch_[j]);
			worker_->setColBounds(branch_row_to_col_[nrows_orig_ + j], rlbd_branch_[j], rubd_branch_[j]);
		}
		/** apply new bounds */
		for (unsigned j = 0; j < branchobj->index_.size(); ++j) {
			si_->setRowBounds(branchobj->index_[j], branchobj->lb_[j], branchobj->ub_[j]);
			worker_->setColBounds(branch_row_to_col_[branchobj->index_[j]], branchobj->lb_[j], branchobj->ub_[j]);
			DSPdebugMessage("Set subproblem column bounds: %d [%e %e]\n",
					branch_row_to_col_[branchobj->index_[j]], branchobj->lb_[j], branchobj->ub_[j]);
		}
	}
}
