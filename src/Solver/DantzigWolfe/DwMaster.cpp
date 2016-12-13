//
// Created by Kibaek Kim on 8/27/16.
//

//#define DSP_DEBUG

#include "cplex.h"
/** Coin */
#include "OsiClpSolverInterface.hpp"
#include "OsiCpxSolverInterface.hpp"
#include "CoinUtility.hpp"
/** Dsp */
#include "Utility/DspUtility.h"
#include "Solver/DantzigWolfe/DwMaster.h"
#include "Model/TssModel.h"

DwMaster::DwMaster(DwWorker* worker):
DecSolver(worker->model_, worker->par_, worker->message_),
useCpxBarrier_(false),
phase_(1),
rlbd_branch_(NULL),
rubd_branch_(NULL),
worker_(worker),
ncols_orig_(0),
ncols_start_(0),
nrows_(0),
nrows_orig_(0),
nrows_branch_(0),
nrows_conv_(0),
org_mat_(NULL),
org_clbd_(NULL),
org_cubd_(NULL),
org_obj_(NULL),
org_ctype_(NULL),
org_rlbd_(NULL),
org_rubd_(NULL),
itercnt_(0) {}

DwMaster::~DwMaster() {
	for (unsigned i = 0; i < auxcols_.size(); ++i)
		FREE_PTR(auxcols_[i]);
	FREE_ARRAY_PTR(rlbd_branch_);
	FREE_ARRAY_PTR(rubd_branch_);
	FREE_PTR(org_mat_);
	FREE_ARRAY_PTR(org_clbd_);
	FREE_ARRAY_PTR(org_cubd_);
	FREE_ARRAY_PTR(org_obj_);
	FREE_ARRAY_PTR(org_ctype_);
	FREE_ARRAY_PTR(org_rlbd_);
	FREE_ARRAY_PTR(org_rubd_);
	for (unsigned i = 0; i < cols_generated_.size(); ++i)
		FREE_PTR(cols_generated_[i]);
}

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
	nrows_conv_ = worker_->getNumSubprobs(); /**< number of convex combination rows in the restricted master */

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

	OsiCpxSolverInterface* cpx = NULL;

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
	//si_->passInMessageHandler(si_->messageHandler());
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
	useCpxBarrier_ = par_->getBoolParam("DW/MASTER/IPM");
	cpx = dynamic_cast<OsiCpxSolverInterface*>(si_);

	if (useCpxBarrier_ && cpx) {
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

DSP_RTN_CODE DwMaster::solve() {
	BGN_TRY_CATCH

	dualobj_ = -COIN_DBL_MAX;

	itercnt_ = 0;

	if (phase_ == 2) {
		DSP_RTN_CHECK_RTN_CODE(solvePhase2());
		if (status_ == DSP_STAT_PRIM_INFEASIBLE) {
			DSPdebugMessage("Converting to Phase 1.\n");
			DSP_RTN_CHECK_RTN_CODE(solvePhase1());
		}
	} else
		DSP_RTN_CHECK_RTN_CODE(solvePhase1());

	if (phase_ == 1) {
		if (status_ == DSP_STAT_FEASIBLE || status_ == DSP_STAT_OPTIMAL) {
			if (primobj_ > 1.0e-8)
				status_ = DSP_STAT_PRIM_INFEASIBLE;
			else {
				DSPdebugMessage("Converting to Phase 2.\n");
				DSP_RTN_CHECK_RTN_CODE(solvePhase2());
			}
		}
	}

	/** collect solutions */
	/** collect solutions */
	if (status_ == DSP_STAT_OPTIMAL ||
		status_ == DSP_STAT_FEASIBLE ||
		status_ == DSP_STAT_LIM_ITERorTIME) {
		/** recover original solution */
		CoinZeroN(primsol_, ncols_orig_);
		for (unsigned k = 0, j = 0; k < cols_generated_.size(); ++k) {
			/** do not consider inactive columns */
			if (cols_generated_[k]->active_ == false)
				continue;
			CoinPackedVector xlam = cols_generated_[k]->x_ * si_->getColSolution()[j];
			for (int i = 0; i < xlam.getNumElements(); ++i)
				primsol_[xlam.getIndices()[i]] += xlam.getElements()[i];
			j++;
		}

		/** heuristics */
		DSP_RTN_CHECK_RTN_CODE(heuristics());
	} else if (status_ == DSP_STAT_DUAL_INFEASIBLE) {
		//DSPdebug(si_->writeMps("master"));
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DwMaster::initialColumns() {
#define FREE_MEMORY     \
	FREE_ARRAY_PTR(piA) \
	FREE_ARRAY_PTR(Ax)

	/** column generation info */
	double* piA = NULL;
	double* Ax = NULL;
	std::vector<int> subinds;
	std::vector<int> substatus;
	std::vector<double> subcxs;
	std::vector<double> subobjs;
	std::vector<CoinPackedVector*> subsols;

	BGN_TRY_CATCH

	/** allocate memory */
	piA = new double [ncols_orig_];
	Ax  = new double [nrows_orig_];

	/** generate columns */
	CoinZeroN(piA, ncols_orig_);
	worker_->generateCols(2, piA, subinds, substatus, subcxs, subobjs, subsols);
	dualobj_ = 0.0;

	for (unsigned int s = 0; s < subinds.size(); ++s) {

		int sind = subinds[s];                  /**< actual subproblem index */
		const CoinPackedVector* x = subsols[s]; /**< retrieve subproblem solution */
		double newcoef = subobjs[s];            /**< column objective coefficient */
		//DSPdebugMessage("subproblem %d objective value %e\n", sind, newcoef);

		/** calculate initial dual bound */
		dualobj_ += subobjs[s];

		/** take A x^k */
		org_mat_->times(*x, Ax);

		/** create a column vector */
		CoinPackedVector colvec;
		colvec.reserve(nrows_);

		/** original constraints */
		for (int i = 0; i < nrows_orig_; ++i) {
			if (fabs(Ax[i]) > 1.0e-10)
				colvec.insert(i, Ax[i]);
		}

		/** branching constraints */
		for (int i = 0; i < nrows_branch_; ++i) {
			int j = branch_row_to_col_[nrows_orig_+i];
			int sparse_index = x->findIndex(j);
			if (sparse_index == -1) continue;
			double val = x->getElements()[sparse_index];
			if (fabs(val) > 1.0e-10)
				colvec.insert(nrows_orig_+i, val);
		}

		/** convex combination constraints, unless unbounded */
		if (substatus[s] != DSP_STAT_DUAL_INFEASIBLE)
			colvec.insert(nrows_orig_ + nrows_branch_ + sind, 1.0);
		//DSPdebugMessage("new column from subproblem %d:\n", sind);
		//DSPdebug(DspMessage::printArray(&colvec));

		/** store columns */
		CoinPackedVector copy_of_x(*x);
		cols_generated_.push_back(new DwCol(sind, copy_of_x, colvec, newcoef, 0.0, COIN_DBL_MAX));
	}
	message_->print(0, "Generated %u initial columns. Initial dual bound %e\n", cols_generated_.size(), dualobj_);

	/** clear column generation information */
	for (unsigned i = 0; i < subsols.size(); ++i)
		FREE_PTR(subsols[i]);

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

DSP_RTN_CODE DwMaster::solvePhase1() {
	BGN_TRY_CATCH

	if (phase_ == 2) {
		/** set phase 1 problem */
		int ncols = si_->getNumCols();
		for (int j = 0; j < ncols; ++j)
			si_->setObjCoeff(j, 0.0);
		for (unsigned j = 0; j < auxcols_.size(); ++j) {
			si_->addCol(*auxcols_[j], 0.0, COIN_DBL_MAX, 1.0);
			auxcolindices_.push_back(ncols + j);
		}
		phase_ = 1;
		DSPdebugMessage("Phase 1 has %d rows and %d columns.\n", si_->getNumRows(), si_->getNumCols());
	}

	DSP_RTN_CHECK_RTN_CODE(gutsOfSolve());

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

DSP_RTN_CODE DwMaster::solvePhase2() {
	BGN_TRY_CATCH

	if (phase_ == 1) {
		/** set phase 2 problem */
		//DSPdebug(si_->writeMps("beforePhase2"));
		si_->deleteCols(auxcolindices_.size(), &auxcolindices_[0]);
		auxcolindices_.clear();
		DSPdebugMessage("Phase 2 has %d rows and %d columns.\n", si_->getNumRows(), si_->getNumCols());
		for (unsigned k = 0, j = 0; k < cols_generated_.size(); ++k)
			if (cols_generated_[k]->active_) {
				if (j >= si_->getNumCols()) {
					message_->print(0, "Trying to access invalid column index %d (ncols %d)\n", j, si_->getNumCols());
					return DSP_RTN_ERR;
				}
				si_->setObjCoeff(j, cols_generated_[k]->obj_);
				j++;
			}
		phase_ = 2;
		//DSPdebug(si_->writeMps("initialPhase2"));
	}

	DSP_RTN_CHECK_RTN_CODE(gutsOfSolve());
	//DSPdebug(si_->writeMps("afterPhase2"));

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

DSP_RTN_CODE DwMaster::gutsOfSolve() {
#define FREE_MEMORY       \
	FREE_PTR(ws)          \
	FREE_ARRAY_PTR(price) \
	FREE_ARRAY_PTR(piA)

	OsiCpxSolverInterface* cpx = NULL;
	CoinWarmStartBasis* ws = NULL;
	double* price = NULL;
	double* piA = NULL;
	bool hasBasis = false;

	double feastol = 1.0e-8; /**< feasibility tolerance */
	double curLb = -COIN_DBL_MAX;
	double relgap;

	/** column generation info */
	int nColsAdded = 0;

	BGN_TRY_CATCH

	cpx = dynamic_cast<OsiCpxSolverInterface*>(si_);

	/** allocate memory */
	price = new double [si_->getNumRows()];
	piA = new double [ncols_orig_];

	std::vector<double> prevsol(si_->getNumCols(), 0.0);

	/** timing results */
	double stime;
	double t_total = CoinGetTimeOfDay();
	double t_master = 0.0;
	double t_gencols = 0.0;

	/** use dual simplex after branching */
	if (!useCpxBarrier_ || cpx == NULL) {
		si_->setHintParam(OsiDoDualInResolve, true);
		si_->messageHandler()->setLogLevel(0);
		hasBasis = true;
	}

	/** resolve */
	stime = CoinGetTimeOfDay();
	si_->resolve();
	t_master += CoinGetTimeOfDay() - stime;
	convertCoinToDspStatus(si_, status_);

	/** calculate primal objective value */
	DSP_RTN_CHECK_RTN_CODE(calculatePrimalObjective());

	/** get price */
	CoinCopyN(si_->getRowPrice(), si_->getNumRows(), price);

	/** get warm-start information */
	if (hasBasis)
		DSP_RTN_CHECK_RTN_CODE(getWarmStartInfo(prevsol, ws));

	/** use primal simplex after column generation */
	if (!useCpxBarrier_)
		si_->setHintParam(OsiDoDualInResolve, false);

	while (status_ == DSP_STAT_OPTIMAL) {

		/** column management */
		DSP_RTN_CHECK_RTN_CODE(reduceCols());

		/** generate columns */
		stime = CoinGetTimeOfDay();
		DSP_RTN_CHECK_RTN_CODE(generateCols(price, piA, curLb, nColsAdded));
		t_gencols += CoinGetTimeOfDay() - stime;

		/** subproblem solution may declare infeasibility. */
		if (status_ == DSP_STAT_PRIM_INFEASIBLE)
			break;

		/** update master */
		DSP_RTN_CHECK_RTN_CODE(updateModel(price, curLb));

		/** set warm-start information */
		if (hasBasis)
			DSP_RTN_CHECK_RTN_CODE(setWarmStartInfo(prevsol, ws));

		/** re-optimize the master */
		stime = CoinGetTimeOfDay();
		si_->resolve();
		t_master += CoinGetTimeOfDay() - stime;
		convertCoinToDspStatus(si_, status_);

		/** calculate primal objective value */
		DSP_RTN_CHECK_RTN_CODE(calculatePrimalObjective());

		/** get price */
		CoinCopyN(si_->getRowPrice(), si_->getNumRows(), price);
#ifdef DSP_DEBUG
		double tr_terms = 0.0;
		for (int j = 0; j < 2 * (nrows_branch_ + nrows_orig_); ++j)
			tr_terms += si_->getObjCoefficients()[j] * si_->getColSolution()[j];
		DSPdebugMessage("master obj %e, obj w/o tr %e\n", primobj_, primobj_ - tr_terms);
#endif

		relgap = (primobj_-dualobj_)/(1.0e-10+fabs(primobj_))*100;
		message_->print(2, "[Phase %d] Iteration %3d: Master objective %e, ", phase_, itercnt_, primobj_);
		if (phase_ == 2)
			message_->print(2, "Lb %e (gap %.2f %%), ", dualobj_, relgap);
		message_->print(2, "nrows %d, ncols %d (new cols %d), ", si_->getNumRows(), si_->getNumCols(), nColsAdded);
		if (!useCpxBarrier_)
			message_->print(2, "itercnt %d, ", si_->getIterationCount());
		message_->print(2, "timing (total %.2f, master %.2f, gencols %.2f), statue %d\n",
				CoinGetTimeOfDay() - t_total, t_master, t_gencols, status_);

		/** get warm-start information */
		if (hasBasis)
			DSP_RTN_CHECK_RTN_CODE(getWarmStartInfo(prevsol, ws));

		/** termination test */
		if (terminationTest(nColsAdded, itercnt_, relgap))
			break;

		itercnt_++;

		/** run heuristics */
//		DSP_RTN_CHECK_RTN_CODE(heuristics());

#ifdef DSP_DEBUG1
		char fname[128];
		sprintf(fname, "master%d", itercnt);
		si_->writeMps(fname);
#endif
	}

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

DSP_RTN_CODE DwMaster::calculatePrimalObjective() {
	BGN_TRY_CATCH
	primobj_ = si_->getObjValue();
	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

DSP_RTN_CODE DwMaster::restoreCols() {
	BGN_TRY_CATCH
	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

DSP_RTN_CODE DwMaster::reduceCols() {
	BGN_TRY_CATCH
	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

DSP_RTN_CODE DwMaster::generateCols(
		const double* price, /**< [in] price */
		double*& piA,        /**< [out] pi^T A */
		double& lb,          /**< [out] lower bound (only for phase 2) */
		int& ncols           /**< [out] number of columns generated */) {
	/** column generation info */
	std::vector<int> subinds;
	std::vector<int> substatus;
	std::vector<double> subcxs;
	std::vector<double> subobjs;
	std::vector<CoinPackedVector*> subsols;

	BGN_TRY_CATCH

	/** calculate pi^T A */
	DSP_RTN_CHECK_RTN_CODE(calculatePiA(price, piA));

	/** generate columns */
	DSP_RTN_CHECK_RTN_CODE(
			worker_->generateCols(phase_, piA, subinds, substatus, subcxs, subobjs, subsols));

	/** termination test */
	terminationTestColgen(substatus);

	if (status_ == DSP_STAT_FEASIBLE) {
		/** calculate lower bound */
		if (phase_ == 2) {
			DSP_RTN_CHECK_RTN_CODE(getLagrangianBound(price, subobjs, lb));
			DSPdebugMessage("Current lower bound %e, best lower bound %e\n", lb, dualobj_);
		}

		/** create and add columns */
		DSP_RTN_CHECK_RTN_CODE(
				addCols(price, piA, subinds, substatus, subcxs, subobjs, subsols, ncols));
	}

	/** free memory for subproblem solutions */
	for (unsigned i = 0; i < subsols.size(); ++i)
		FREE_PTR(subsols[i]);
	subsols.clear();

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

DSP_RTN_CODE DwMaster::addCols(
		const double* price,                  /**< [in] price */
		const double* piA,                    /**< [in] pi^T A */
		std::vector<int>& indices,            /**< [in] subproblem indices corresponding to cols*/
		std::vector<int>& statuses,           /**< [in] subproblem solution status */
		std::vector<double>& cxs,             /**< [in] solution times original objective coefficients */
		std::vector<double>& objs,            /**< [in] subproblem objective values */
		std::vector<CoinPackedVector*>& sols, /**< [in] subproblem solutions */
		int& nadded                           /**< [out] number of columns added */) {
#define FREE_MEMORY \
	FREE_ARRAY_PTR(Ax)

	double* Ax = NULL;
	CoinPackedVector colvec;
	TssModel* tss = NULL;

	BGN_TRY_CATCH

	/** allocate memory */
	Ax = new double [nrows_orig_];
	colvec.reserve(nrows_);

	if (model_->isStochastic())
		tss = dynamic_cast<TssModel*>(model_);

	/** reset counter */
	nadded = 0;

	for (unsigned int s = 0; s < indices.size(); ++s) {
		int sind = indices[s]; /**< actual subproblem index */

		/** cutoff = dual variable corresponding to the convex-combination constraint */
		double cutoff = price[nrows_orig_ + nrows_branch_ + sind];
		DSPdebugMessage("pricing out: %e < %e ? (colobj %e, status %d)\n", objs[s], cutoff, cxs[s], statuses[s]);

		if (statuses[s] == DSP_STAT_DUAL_INFEASIBLE || objs[s] < cutoff - 1.0e-4) {
			/** retrieve subproblem solution */
			const CoinPackedVector* x = sols[s];

			/** create a column objective */
			double newcoef = cxs[s];

			/** take A x^k */
			org_mat_->times(*x, Ax);

			/** clear a column vector */
			colvec.clear();

			/** original constraints */
			for (int i = 0; i < nrows_orig_; ++i)
				if (fabs(Ax[i]) > 1.0e-10)
					colvec.insert(i, Ax[i]);

			/** skip zero column */
//			if (colvec.getNumElements() == 0)
//				continue;

			/** branching constraints */
			for (int i = 0; i < nrows_branch_; ++i) {
				int j = branch_row_to_col_[nrows_orig_+i];
				int sparse_index = x->findIndex(j);
				if (sparse_index == -1) continue;
				double val = x->getElements()[sparse_index];
				if (fabs(val) > 1.0e-10)
					colvec.insert(nrows_orig_ + i, val);
			}

			/** convex combination constraints */
			if (statuses[s] != DSP_STAT_DUAL_INFEASIBLE)
				colvec.insert(nrows_orig_ + nrows_branch_ + sind, 1.0);
//			DSPdebugMessage("new column from subproblem %d (objcoef %e):\n", sind, newcoef);
//			DSPdebug(DspMessage::printArray(&colvec));

			/** add the column vector */
			if (phase_ == 1)
				si_->addCol(colvec, 0.0, COIN_DBL_MAX, 0.0);
			else if (phase_ == 2)
				si_->addCol(colvec, 0.0, COIN_DBL_MAX, newcoef);

			/** store columns */
			cols_generated_.push_back(new DwCol(sind, *x, colvec, newcoef, 0.0, COIN_DBL_MAX));
			nadded++;
		}
	}
	//message_->print(3, "Number of columns in the pool: %u\n", cols_generated_.size());

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

DSP_RTN_CODE DwMaster::getLagrangianBound(
		const double* price,       /**< [in] price */
		std::vector<double>& objs, /**< [in] subproblem objective values */
		double& lb                 /**< [out] lower bound */) {
	BGN_TRY_CATCH

	/** calculate lower bound */
#if 1
	lb = 0.0;
	const double* rlbd = si_->getRowLower();
	const double* rubd = si_->getRowUpper();
	for (int j = 0; j < nrows_orig_ + nrows_branch_; ++j) {
		if (fabs(price[j]) < 1.0e-10)
			continue;
		else if (price[j] > 0)
			lb += price[j] * rlbd[j];
		else
			lb += price[j] * rubd[j];
	}
	DSPdebugMessage("pia %e\n", lb);
#else
	//lb = primobj_;
	lb = si_->getObjValue();
	for (int j = nrows_orig_ + nrows_branch_; j < nrows_; ++j)
		lb -= price[j];
	DSPdebugMessage("si_->getObjValue() %e price %e\n", si_->getObjValue(), lb - primobj_);
#endif
	for (unsigned int s = 0; s < objs.size(); ++s) {
		lb += objs[s];
		DSPdebugMessage("subobj[%d] %e\n", s, objs[s]);
	}
	DSPdebugMessage("lb %e\n", lb);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

bool DwMaster::terminationTest(int nnewcols, int itercnt, double relgap) {

	if (nnewcols == 0)
		return true;

	bool term = false;
	double feastol = 1.0e-8;

	if (phase_ == 1) {
		if (si_->isProvenOptimal() && primobj_ < feastol) {
			status_ = DSP_STAT_FEASIBLE;
			DSPdebugMessage("Phase 1 found a feasible solution! %e < %e\n", si_->getObjValue(), feastol);
			term = true;
		}
	}
	if (phase_ == 2) {
		if (iterlim_ <= itercnt /*||
			relgap < 0.01 ||
			dualobj_ >= bestprimobj_*/) {
			status_ = DSP_STAT_FEASIBLE;
			term = true;
		}
	}

	return term;
}

DSP_RTN_CODE DwMaster::updateModel(
		const double* price, /**< [in] price */
		double curLb         /**< [in] current lower bound */) {
	BGN_TRY_CATCH

	/** update the best dual objective */
	if (dualobj_ < curLb) dualobj_ = curLb;

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

DSP_RTN_CODE DwMaster::calculatePiA(
		const double* price, /**< [in] price */
		double*& piA         /**< [out] pi^T A */) {
	BGN_TRY_CATCH

	org_mat_->transposeTimes(price, piA);
	for (int i = nrows_orig_; i < nrows_orig_ + nrows_branch_; ++i)
		piA[branch_row_to_col_[i]] += price[i];
//	DSPdebugMessage("piA (%d):\n", ncols_orig_);
//	DSPdebug(DspMessage::printArray(ncols_orig_, piA));

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

bool DwMaster::terminationTestColgen(std::vector<int>& statuses) {
	bool term = false;
	status_ = DSP_STAT_FEASIBLE;
	for (unsigned s = 0; s < statuses.size(); ++s) {
		if (statuses[s] != DSP_STAT_OPTIMAL && statuses[s] != DSP_STAT_DUAL_INFEASIBLE) {
			status_ = DSP_STAT_PRIM_INFEASIBLE;
			break;
		}
		DSPdebugMessage("statuses[%d] %d\n", s, statuses[s]);
	}
	if (status_ == DSP_STAT_PRIM_INFEASIBLE) {
		DSPdebugMessage("Subproblems are infeasible.\n");
		term = true;
	}
	return term;
}

DSP_RTN_CODE DwMaster::getWarmStartInfo(
		std::vector<double>& sol, /**< [out] current solution */
		CoinWarmStartBasis*& ws   /**< [out] warmstart basis */) {
	BGN_TRY_CATCH

	/** store previous solution and basis */
	CoinCopyN(si_->getColSolution(), si_->getNumCols(), &sol[0]);
	ws = dynamic_cast<CoinWarmStartBasis*>(si_->getWarmStart());

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DwMaster::setWarmStartInfo(
		std::vector<double>& sol, /**< [out] current solution */
		CoinWarmStartBasis*& ws   /**< [out] warmstart basis */) {
	BGN_TRY_CATCH

	/** set previous solution and basis */
	sol.resize(si_->getNumCols(), 0.0);
	si_->setColSolution(&sol[0]);

	ws->resize(si_->getNumRows(), si_->getNumCols());
	si_->setWarmStart(ws);
	FREE_PTR(ws);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DwMaster::heuristics() {
	if (!runHeuristics_)
		return DSP_RTN_OK;

	double stime;

	BGN_TRY_CATCH

	if (iterlim_ == 1 && par_->getBoolParam("DW/HEURISTICS/TRIVIAL")) {
		message_->print(1, "Heuristic (trivial) searches solutions...\n");
		stime = CoinGetTimeOfDay();
		DSP_RTN_CHECK_RTN_CODE(heuristicTrivial());
		message_->print(1, "Heuristic (trivial) spent %.2f seconds [best %e].\n", CoinGetTimeOfDay() - stime, bestprimobj_);
	}

	if (par_->getBoolParam("DW/HEURISTICS/FP1")) {
		message_->print(1, "Heuristic (FP-like[+1]) searches solutions....\n");
		stime = CoinGetTimeOfDay();
		DSP_RTN_CHECK_RTN_CODE(heuristicFp(1));
		message_->print(1, "Heuristic (FP-like[+1]) spent %.2f seconds [best %e].\n", CoinGetTimeOfDay() - stime, bestprimobj_);
	}

	if (par_->getBoolParam("DW/HEURISTICS/FP2")) {
		message_->print(1, "Heuristic (FP-like[-1]) searches solutions....\n");
		stime = CoinGetTimeOfDay();
		DSP_RTN_CHECK_RTN_CODE(heuristicFp(-1));
		message_->print(1, "Heuristic (FP-like[-1]) spent %.2f seconds [best %e].\n", CoinGetTimeOfDay() - stime, bestprimobj_);
	}

	if (par_->getBoolParam("DW/HEURISTICS/DIVE")) {
		message_->print(1, "Heuristic (dive) searches solutions...\n");
		stime = CoinGetTimeOfDay();
		DSP_RTN_CHECK_RTN_CODE(heuristicDive());
		message_->print(1, "Heuristic (dive) spent %.2f seconds [best %e].\n", CoinGetTimeOfDay() - stime, bestprimobj_);
	}

	/** restore the original settings */
	iterlim_ = par_->getIntParam("DW/ITER_LIM");

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

DSP_RTN_CODE DwMaster::preHeuristic(
		double*& rlbd,    /**< [out] original row lower bounds */
		double*& rubd,    /**< [out] original row lower bounds */
		double*& primsol, /**< [out] original primal solution */
		double& primobj,  /**< [out] original primal objective */
		double& dualobj,  /**< [out] original dual objective */
		int& status       /**< [out] original solution status */) {
	BGN_TRY_CATCH

#ifndef DSP_DEBUG
	message_->logLevel_ = 1;
#endif

	/** save the original row bounds */
	rlbd = new double [nrows_branch_];
	rubd = new double [nrows_branch_];
	CoinCopyN(si_->getRowLower() + nrows_orig_, nrows_branch_, rlbd);
	CoinCopyN(si_->getRowUpper() + nrows_orig_, nrows_branch_, rubd);

	/** save the original solutions */
	primsol = new double [ncols_orig_];
	CoinCopyN(primsol_, ncols_orig_, primsol);
	primobj = primobj_;
	dualobj = dualobj_;
	status = status_;

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DwMaster::postHeuristic(
		double*& rlbd,    /**< [out] original row lower bounds */
		double*& rubd,    /**< [out] original row lower bounds */
		double*& primsol, /**< [out] original primal solution */
		double& primobj,  /**< [out] original primal objective */
		double& dualobj,  /**< [out] original dual objective */
		int& status       /**< [out] original solution status */) {
	BGN_TRY_CATCH

	message_->logLevel_ = par_->getIntParam("LOG_LEVEL");

	/** restore the original solutions */
	CoinCopyN(primsol, ncols_orig_, primsol_);
	primobj_ = primobj;
	dualobj_ = dualobj;
	status_ = status;

	/** restore the original row bounds */
	std::vector<int> branchIndices;
	for (int i = 0; i < nrows_branch_; ++i) {
		si_->setRowBounds(nrows_orig_ + i, rlbd[i], rubd[i]);
		branchIndices.push_back(branch_row_to_col_[nrows_orig_ + i]);
	}
	worker_->setColBounds(nrows_branch_, &branchIndices[0], rlbd, rubd);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DwMaster::heuristicTrivial() {
#define FREE_MEMORY \
	FREE_ARRAY_PTR(rlbd) \
	FREE_ARRAY_PTR(rubd) \
	FREE_ARRAY_PTR(primsol_org)

	/** row bounds for branching rows */
	double* rlbd = NULL;
	double* rubd = NULL;

	/** original solutions */
	double* primsol_org = NULL;
	double primobj_org, dualobj_org;
	int status_org;
	int iterlim;

	BGN_TRY_CATCH

	DSP_RTN_CHECK_RTN_CODE(
			preHeuristic(rlbd, rubd, primsol_org, primobj_org, dualobj_org, status_org));

	/** fix bounds */
	std::vector<int> branchIndices;
	std::vector<double> branchBounds;
	for (int i = 0; i < nrows_branch_; ++i) {
		/** do not consider those with negative objective coefficients */
		if (org_obj_[branch_row_to_col_[nrows_orig_ + i]] >= 0) continue;
		/** skip fixed bounds */
		if (rlbd[i] == rubd[i]) continue;
		si_->setRowBounds(nrows_orig_ + i, 0.0, 0.0);
		branchIndices.push_back(branch_row_to_col_[nrows_orig_ + i]);
		branchBounds.push_back(0.0);
	}
	worker_->setColBounds(branchIndices.size(), &branchIndices[0], &branchBounds[0], &branchBounds[0]);

	/** solve */
	itercnt_ = 0;
	iterlim_ = par_->getIntParam("DW/HEURISTICS/TRIVIAL/ITER_LIM");
	time_remains_ = par_->getDblParam("DW/HEURISTICS/TRIVIAL/TIME_LIM");
	DSP_RTN_CHECK_RTN_CODE(solvePhase1());

	if (primobj_ > 1.0e-8)
		status_ = DSP_STAT_PRIM_INFEASIBLE;
	else {
		itercnt_ = 0;
		DSP_RTN_CHECK_RTN_CODE(solvePhase2());
	}

	/** collect solutions */
	if (status_ == DSP_STAT_OPTIMAL ||
			status_ == DSP_STAT_FEASIBLE ||
			status_ == DSP_STAT_LIM_ITERorTIME) {
		/** check integer feasibility */
		bool fractional = false;
		for (int i = 0; i < nrows_branch_; ++i) {
			double x = si_->getRowActivity()[nrows_orig_ + i];
			if (fabs(x - floor(x + 0.5)) > 1.0e-10) {
				fractional = true;
				DSPdebugMessage("Heuristic found a fractional solution.\n");
				break;
			}
		}
		if (!fractional) {
			message_->print(1, "Heuristic found an integer solution (objective %e).\n", primobj_);
			if (bestprimobj_ > primobj_) {
				bestprimobj_ = primobj_;
				/** recover original solution */
				CoinZeroN(bestprimsol_, ncols_orig_);
				for (unsigned k = 0, j = ncols_start_; k < cols_generated_.size(); ++k) {
					/** do not consider inactive columns */
					if (cols_generated_[k]->active_ == false)
						continue;
					CoinPackedVector xlam = cols_generated_[k]->x_ * si_->getColSolution()[j];
					for (int i = 0; i < xlam.getNumElements(); ++i)
						bestprimsol_[xlam.getIndices()[i]] += xlam.getElements()[i];
					j++;
				}
				message_->print(1, "Heuristic updated the best upper bound %e.\n", bestprimobj_);
			}
		}
	} else if (status_ == DSP_STAT_DUAL_INFEASIBLE) {
		//DSPdebug(si_->writeMps("master"));
	}

	DSP_RTN_CHECK_RTN_CODE(
			postHeuristic(rlbd, rubd, primsol_org, primobj_org, dualobj_org, status_org));

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

DSP_RTN_CODE DwMaster::heuristicFp(int direction) {
#define FREE_MEMORY \
	FREE_ARRAY_PTR(rlbd) \
	FREE_ARRAY_PTR(rubd) \
	FREE_ARRAY_PTR(primsol_org)

	/** row bounds for branching rows */
	double* rlbd = NULL;
	double* rubd = NULL;

	/** original solutions */
	double* primsol_org = NULL;
	double primobj_org, dualobj_org;
	int status_org;

	BGN_TRY_CATCH

	DSP_RTN_CHECK_RTN_CODE(
			preHeuristic(rlbd, rubd, primsol_org, primobj_org, dualobj_org, status_org));

	/** round and fix bounds */
	for (int i = 0; i < nrows_branch_; ++i) {
		/** do not consider those with non-positive objective coefficients */
		int colind = branch_row_to_col_[nrows_orig_ + i];
		if (org_obj_[colind] * direction <= 0) continue;

		double rounded = round(primsol_[colind]);
		si_->setRowBounds(nrows_orig_ + i, rounded, rounded);
		worker_->setColBounds(1, &colind, &rounded, &rounded);
	}

	/** solve */
	DSP_RTN_CHECK_RTN_CODE(solvePhase1());

	if (primobj_ > 1.0e-8)
		status_ = DSP_STAT_PRIM_INFEASIBLE;
	else
		DSP_RTN_CHECK_RTN_CODE(solvePhase2());

	/** collect solutions */
	if (status_ == DSP_STAT_OPTIMAL || status_ == DSP_STAT_FEASIBLE) {
		/** check integer feasibility */
		bool fractional = false;
		for (int i = 0; i < nrows_branch_; ++i) {
			double x = si_->getRowActivity()[nrows_orig_ + i];
			if (fabs(x - floor(x + 0.5)) > 1.0e-10) {
				fractional = true;
				DSPdebugMessage("Heuristic found a fractional solution.\n");
				break;
			}
		}
		if (!fractional) {
			DSPdebugMessage("Heuristic found an integer solution.\n");
			if (bestprimobj_ > primobj_) {
				bestprimobj_ = primobj_;
				/** recover original solution */
				CoinZeroN(bestprimsol_, ncols_orig_);
				for (unsigned k = 0, j = ncols_start_; k < cols_generated_.size(); ++k) {
					/** do not consider inactive columns */
					if (cols_generated_[k]->active_ == false)
						continue;
					CoinPackedVector xlam = cols_generated_[k]->x_ * si_->getColSolution()[j];
					for (int i = 0; i < xlam.getNumElements(); ++i)
						bestprimsol_[xlam.getIndices()[i]] += xlam.getElements()[i];
					j++;
				}
				DSPdebugMessage("Heuristic updated the best upper bound %e.\n", bestprimobj_);
			}
		}
	} else if (status_ == DSP_STAT_DUAL_INFEASIBLE) {
		//DSPdebug(si_->writeMps("master"));
	}

	DSP_RTN_CHECK_RTN_CODE(
			postHeuristic(rlbd, rubd, primsol_org, primobj_org, dualobj_org, status_org));

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

DSP_RTN_CODE DwMaster::heuristicDive() {
#define FREE_MEMORY \
	FREE_ARRAY_PTR(rlbd) \
	FREE_ARRAY_PTR(rubd) \
	FREE_ARRAY_PTR(primsol_org)

	/** row bounds for branching rows */
	double* rlbd = NULL;
	double* rubd = NULL;

	/** original solutions */
	double* primsol_org = NULL;
	double primobj_org, dualobj_org;
	int status_org;

	std::vector<CoinTriple<int,int,double> > branchList;

	BGN_TRY_CATCH

	DSP_RTN_CHECK_RTN_CODE(
			preHeuristic(rlbd, rubd, primsol_org, primobj_org, dualobj_org, status_org));

	DSP_RTN_CHECK_RTN_CODE(gutsOfDive(branchList, 0));

	DSP_RTN_CHECK_RTN_CODE(
			postHeuristic(rlbd, rubd, primsol_org, primobj_org, dualobj_org, status_org));

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

DSP_RTN_CODE DwMaster::gutsOfDive(
		std::vector<CoinTriple<int, int, double> > branchList, int depth) {

	/** parameters for dive with backtracking */
	int maxdiscrepancy = 2;
	int maxdepth = 3;

	int findPhase = 0;
	int branchIndex, branchDirection, solvePhase;
	double branchValue;
	double dist, mindist;
	bool isInteger;

	int status; /**< status at the current depth */

	BGN_TRY_CATCH

	while (1) {

		if (depth > 0) {
			/** solve */
			itercnt_ = 0;
			iterlim_ = par_->getIntParam("DW/HEURISTICS/DIVE/ITER_LIM");
			time_remains_ = par_->getDblParam("DW/HEURISTICS/DIVE/TIME_LIM");
			DSP_RTN_CHECK_RTN_CODE(solvePhase1());

			/** recover original solution */
			CoinZeroN(primsol_, ncols_orig_);
			for (unsigned k = 0, j = ncols_start_; k < cols_generated_.size(); ++k) {
				/** do not consider inactive columns */
				if (cols_generated_[k]->active_ == false)
					continue;
				CoinPackedVector xlam = cols_generated_[k]->x_ * si_->getColSolution()[j];
				for (int i = 0; i < xlam.getNumElements(); ++i)
					primsol_[xlam.getIndices()[i]] += xlam.getElements()[i];
				j++;
			}
		}

		/** find a fractional */
		findPhase = 0;
		branchIndex = -1;
		mindist = 1.0;
		isInteger = true;
		while (findPhase < 2 && branchIndex < 0) {
			findPhase++;
			for (int i = 0; i < nrows_branch_; ++i) {
				int colind = branch_row_to_col_[nrows_orig_ + i];
				/** do not consider those with non-negative objective coefficients */
				if (findPhase == 0 && org_obj_[colind] >= 0) continue;
				/** do not consider those with negative objective coefficients */
				if (findPhase == 1 && org_obj_[colind] < 0) continue;

				/** skip if already fixed */
				if (si_->getRowUpper()[nrows_orig_ + i] - si_->getRowLower()[nrows_orig_ + i] < 1.0e-8)
					continue;

				dist = fabs(primsol_[colind] - floor(primsol_[colind] + 0.5));
				if (dist < 1.0e-8) continue;

				/** mark as not integer */
				isInteger = false;

				double candValue = primsol_[colind];
				double candRounded = round(candValue);
				int candDirection = (candRounded < candValue) ? -1 : 1;

				/** check if the branch found is in branchList. */
				for (unsigned j = 0; j < branchList.size(); ++j) {
					if (branchList[j].first == i &&
						branchList[j].second == candDirection &&
						fabs(branchList[j].third - candRounded) < 1.0e-10) {
#if 0
						/** yes, it exists; so do not branch on this. */
						dist = 1.0;
						break;
#else
						/** flip */
						if (candValue > candRounded)
							candRounded += 1.0;
						else
							candRounded -= 1.0;
						dist = 0.5 - dist;
						if (fabs(branchList[j].third - candRounded) < 1.0e-10) {
							dist = 1.0;
							break;
						}
#endif
					}
				}

				if (dist < mindist) {
					mindist = dist;
					branchDirection = candDirection;
					branchIndex = i;
					branchValue = candValue;
				}
			}
		}

		/** found a fractional variable */
		if (branchIndex > -1) {
			/** keep the current node bounds */
			double rlbd_node = si_->getRowLower()[nrows_orig_ + branchIndex];
			double rubd_node = si_->getRowUpper()[nrows_orig_ + branchIndex];

			/** fix bounds */
			double rounded = round(branchValue);
			si_->setRowBounds(nrows_orig_ + branchIndex, rounded, rounded);
			worker_->setColBounds(1, &branch_row_to_col_[nrows_orig_ + branchIndex], &rounded, &rounded);
			message_->print(2, "Diving fixed variable %d [%e] to %e (discrepancy %u, depth %d).\n", branchIndex, branchValue, rounded, branchList.size(), depth);

			/** recursive call */
			status = status_;
			DSP_RTN_CHECK_RTN_CODE(gutsOfDive(branchList, depth+1));
			status_ = status;

			/** restore node bounds */
			si_->setRowBounds(nrows_orig_ + branchIndex, rlbd_node, rubd_node);
			worker_->setColBounds(1, &branch_row_to_col_[nrows_orig_ + branchIndex], &rlbd_node, &rubd_node);

			/** put a branch to the list */
			branchList.push_back(CoinMakeTriple(branchIndex, branchDirection, branchValue));

			DSPdebugMessage("discrepancy %u, depth %d\n", branchList.size(), depth);
			if (branchList.size() > maxdiscrepancy || depth > maxdepth)
				break;
		} else if (!isInteger) {
			break;
		} else {
			message_->print(1, "Diving found an integer solution.\n");
#define FIX_NOW
#ifdef FIX_NOW
			/** backup and fix all bounds */
			std::vector<int> branchIndices;
			std::vector<double> branchBounds;
			double* rlbd_tmp = new double [nrows_branch_];
			double* rubd_tmp = new double [nrows_branch_];
			CoinCopyN(si_->getRowLower() + nrows_orig_, nrows_branch_, rlbd_tmp);
			CoinCopyN(si_->getRowUpper() + nrows_orig_, nrows_branch_, rubd_tmp);
			for (int j = 0; j < nrows_branch_; ++j) {
				double rounded = round(primsol_[branch_row_to_col_[nrows_orig_ + j]]);
				si_->setRowBounds(nrows_orig_ + j, rounded, rounded);
				branchIndices.push_back(branch_row_to_col_[nrows_orig_ + j]);
				branchBounds.push_back(rounded);
			}
			worker_->setColBounds(branchIndices.size(), &branchIndices[0], &branchBounds[0], &branchBounds[0]);

			/** solve */
			itercnt_ = 0;
			iterlim_ = par_->getIntParam("DW/HEURISTICS/DIVE/ITER_LIM");
			time_remains_ = par_->getDblParam("DW/HEURISTICS/DIVE/TIME_LIM");
			DSP_RTN_CHECK_RTN_CODE(solvePhase1());

			/** determine solution status */
			if (primobj_ > 1.0e-8) {
				message_->print(1, "The integer solution is infeasible.\n");
				status_ = DSP_STAT_PRIM_INFEASIBLE;
				break;
			} else
				status_ = DSP_STAT_FEASIBLE;
#endif
			message_->print(1, "Diving is evaluating the integer solution.\n");
			itercnt_ = 0;
			iterlim_ = par_->getIntParam("DW/HEURISTICS/DIVE/ITER_LIM");
			time_remains_ = par_->getDblParam("DW/HEURISTICS/DIVE/TIME_LIM");
			DSP_RTN_CHECK_RTN_CODE(solvePhase2());

			/** collect solutions */
			bool terminateLoop = false;
			if (status_ == DSP_STAT_OPTIMAL || status_ == DSP_STAT_FEASIBLE) {
				/** recover original solution */
				CoinZeroN(primsol_, ncols_orig_);
				for (unsigned k = 0, j = ncols_start_; k < cols_generated_.size(); ++k) {
					/** do not consider inactive columns */
					if (cols_generated_[k]->active_ == false)
						continue;
					CoinPackedVector xlam = cols_generated_[k]->x_ * si_->getColSolution()[j];
					for (int i = 0; i < xlam.getNumElements(); ++i)
						primsol_[xlam.getIndices()[i]] += xlam.getElements()[i];
					j++;
				}
				bool fractional = false;
#ifndef FIX_NOW
				/** check integer feasibility */
				for (int i = 0; i < nrows_branch_; ++i) {
					double x = primsol_[branch_row_to_col_[nrows_orig_ + i]];
					if (fabs(x - floor(x + 0.5)) > 1.0e-10) {
						fractional = true;
						DSPdebugMessage("Heuristic found a fractional solution (x %d [%e]).\n", i, x);
						break;
					}
				}
#endif
				if (!fractional) {
					primobj_ = 0.0;
					for (int j = ncols_start_; j < si_->getNumCols(); ++j)
						primobj_ += si_->getObjCoefficients()[j] * si_->getColSolution()[j];
					message_->print(1, "Diving found an integer solution (objective %e).\n", primobj_);
					if (bestprimobj_ > primobj_) {
						bestprimobj_ = primobj_;
						CoinCopyN(primsol_, ncols_orig_, bestprimsol_);
						message_->print(1, "Diving updated the best upper bound %e.\n", bestprimobj_);
					}
					terminateLoop = true;
				}
			}

#ifdef FIX_NOW
			for (int j = 0; j < nrows_branch_; ++j)
				si_->setRowBounds(nrows_orig_ + j, rlbd_tmp[j], rubd_tmp[j]);
			worker_->setColBounds(branchIndices.size(), &branchIndices[0], rlbd_tmp, rubd_tmp);
			FREE_ARRAY_PTR(rlbd_tmp)
			FREE_ARRAY_PTR(rubd_tmp)
#endif

			if (terminateLoop)
				break;
		}
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
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
	double dist, maxdist = 1.0e-6;
	int branchingIndex = -1;
	double branchingValue;

	BGN_TRY_CATCH

	const double* rlbd = si_->getRowLower();
	const double* rubd = si_->getRowUpper();

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
#define EXPERIMENTAL
		/** most fractional value */
		for (int j = 0; j < ncols_orig_; ++j) {
			if (org_ctype_[j] == 'C') continue;
			if (org_obj_[j] >= 0.0) continue;
			dist = fabs(primsol_[j] - floor(primsol_[j] + 0.5));
			if (dist > maxdist) {
				maxdist = dist;
				branchingIndex = j;
				branchingValue = primsol_[j];
			}
		}

		if (branchingIndex < 0) {
			/** most fractional value */
			for (int j = 0; j < ncols_orig_; ++j) {
				if (org_ctype_[j] == 'C') continue;
				if (org_obj_[j] < 0.0) continue;
				dist = fabs(primsol_[j] - floor(primsol_[j] + 0.5));
				if (dist > maxdist) {
					maxdist = dist;
					branchingIndex = j;
					branchingValue = primsol_[j];
				}
			}
		}
#ifdef EXPERIMENTAL
#else
		/** most fractional value */
		for (int j = 0; j < ncols_orig_; ++j) {
			if (org_ctype_[j] == 'C') continue;
			dist = fabs(primsol_[j] - floor(primsol_[j] + 0.5));
			if (dist > maxdist) {
				maxdist = dist;
				branchingIndex = j;
				branchingValue = primsol_[j];
			}
		}
#endif
	}

	if (branchingIndex > -1) {
		DSPdebugMessage("Creating branch objects on column %d (value %e).\n", branchingIndex, branchingValue);
		branched = true;

		/** creating branching objects */
		branchingUp = new DspBranch();
		branchingDn = new DspBranch();
		for (int j = 0, k = nrows_orig_; j < ncols_orig_; ++j) {
			if (org_ctype_[j] == 'C') continue;
			if (branchingIndex == j) {
				branchingUp->push_back(k, ceil(branchingValue), rubd[k]);
				branchingDn->push_back(k, rlbd[k], floor(branchingValue));
			} else if (rlbd[k] > rlbd_branch_[k - nrows_orig_] || rubd[k] < rubd_branch_[k - nrows_orig_]) {
				/** store any bound changes made in parent nodes */
				branchingUp->push_back(k, rlbd[k], rubd[k]);
				branchingDn->push_back(k, rlbd[k], rubd[k]);
			}
			k++;
		}
	} else {
		DSPdebugMessage("No branch object is found.\n");
		//DSPdebug(si_->writeMps("Incumbent"));
	}

	END_TRY_CATCH_RTN(;,false)

	return branched;
}

void DwMaster::setBranchingObjects(const DspBranch* branchobj) {
#define FREE_MEMORY \
	FREE_ARRAY_PTR(indices)

	int* indices = NULL;

	BGN_TRY_CATCH

	if (branchobj) {
		indices = new int [nrows_branch_];

		/** restore original bounds */
		for (int j = 0; j < nrows_branch_; ++j)
			si_->setRowBounds(nrows_orig_ + j, rlbd_branch_[j], rubd_branch_[j]);
		/** apply new bounds */
		for (unsigned j = 0; j < branchobj->index_.size(); ++j) {
			si_->setRowBounds(branchobj->index_[j], branchobj->lb_[j], branchobj->ub_[j]);
			DSPdebugMessage("Branching on %d [%e, %e]\n",
					branch_row_to_col_[branchobj->index_[j]], branchobj->lb_[j], branchobj->ub_[j]);
		}
		for (int j = 0; j < nrows_branch_; ++j)
			indices[j] = branch_row_to_col_[nrows_orig_ + j];
		worker_->setColBounds(nrows_branch_, indices,
				si_->getRowLower() + nrows_orig_, si_->getRowUpper() + nrows_orig_);
	}

	END_TRY_CATCH(FREE_MEMORY)
	FREE_MEMORY
#undef FREE_MEMORY
}
