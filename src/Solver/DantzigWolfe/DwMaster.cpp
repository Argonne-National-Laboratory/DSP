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
#include "SolverInterface/OoqpEps.h"

DwMaster::DwMaster(DwWorker* worker):
DecSolver(worker->model_, worker->par_, worker->message_),
phase_(1),
worker_(worker),
ncols_orig_(0),
ncols_start_(0),
nrows_(0),
nrows_orig_(0),
nrows_conv_(0),
nrows_core_(0),
nrows_branch_(0),
org_mat_(NULL),
org_clbd_(NULL),
org_cubd_(NULL),
org_obj_(NULL),
org_ctype_(NULL),
org_rlbd_(NULL),
org_rubd_(NULL),
node_clbd_(NULL),
node_cubd_(NULL),
itercnt_(0),
t_total_(0.0),
t_master_(0.0),
t_colgen_(0.0) {
	useBarrier_ = par_->getBoolParam("DW/MASTER/IPM");
}

DwMaster::~DwMaster() {
	FREE_PTR(org_mat_);
	FREE_ARRAY_PTR(org_clbd_);
	FREE_ARRAY_PTR(org_cubd_);
	FREE_ARRAY_PTR(org_obj_);
	FREE_ARRAY_PTR(org_ctype_);
	FREE_ARRAY_PTR(org_rlbd_);
	FREE_ARRAY_PTR(org_rubd_);
	FREE_ARRAY_PTR(node_clbd_);
	FREE_ARRAY_PTR(node_cubd_);
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
	nrows_conv_ = worker_->getNumSubprobs(); /**< number of convex combination rows in the restricted master */
	nrows_core_ = nrows_orig_ + nrows_conv_;
	nrows_branch_ = 0; /**< number of branching rows in the restricted master */

	/** number of rows in the restricted master */
	nrows_ = nrows_core_ + nrows_branch_;

	DSPdebugMessage("nrwos_ %d, nrows_orig_ %d, nrows_conv_ %d, nrows_branch_ %d\n",
			nrows_, nrows_orig_, nrows_conv_, nrows_branch_);

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
	FREE_ARRAY_PTR(rlbd)  \
	FREE_ARRAY_PTR(rubd)

	/** master problem */
	CoinPackedMatrix * mat = NULL;
	double * rlbd = NULL;
	double * rubd = NULL;

	BGN_TRY_CATCH

	/** allocate memory */
	rlbd = new double [nrows_];
	rubd = new double [nrows_];
	node_clbd_ = new double [ncols_orig_];
	node_cubd_ = new double [ncols_orig_];
	CoinCopyN(org_clbd_, ncols_orig_, node_clbd_);
	CoinCopyN(org_cubd_, ncols_orig_, node_cubd_);

	/** create column-wise matrix and set number of rows */
	mat = new CoinPackedMatrix(true, 0, 0);
	mat->setDimensions(nrows_, 0);

	/** Set row bounds */
	CoinCopyN(org_rlbd_, nrows_orig_, rlbd);
	CoinCopyN(org_rubd_, nrows_orig_, rubd);
	CoinFillN(rlbd + nrows_orig_, nrows_conv_, 1.0);
	CoinFillN(rubd + nrows_orig_, nrows_conv_, 1.0);

	/** create solver */
	if (useBarrier_)
		si_ = new OoqpEps();
	else
		si_ = new OsiCpxSolverInterface();
	si_->messageHandler()->logLevel(5);

	/** load problem data */
	si_->loadProblem(*mat, NULL, NULL, NULL, rlbd, rubd);

	DSP_RTN_CHECK_RTN_CODE(initialOsiSolver());

	/** allocate memory for solution */
	bestprimsol_ = new double [ncols_orig_];
	primsol_ = new double [ncols_orig_];

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	/** release memory */
	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

DSP_RTN_CODE DwMaster::initialOsiSolver() {
	BGN_TRY_CATCH
	/** set hint parameters */
	if (useBarrier_) {

	} else {
		OsiCpxSolverInterface* cpx = dynamic_cast<OsiCpxSolverInterface*>(si_);
		if (cpx) {
			CPXsetintparam(cpx->getEnvironmentPtr(), CPX_PARAM_THREADS, 1);
		}
		si_->setHintParam(OsiDoPresolveInResolve, false);
		si_->setHintParam(OsiDoDualInResolve, false);
	}
#if 0
	if (useBarrier_ && cpx) {
		/** use barrier */
		CPXsetintparam(cpx->getEnvironmentPtr(), CPX_PARAM_LPMETHOD, CPX_ALG_BARRIER);
		/** no crossover */
		CPXsetintparam(cpx->getEnvironmentPtr(), CPX_PARAM_BARCROSSALG, -1);
		/** use standard barrier; the others result in numerical instability. */
		CPXsetintparam(cpx->getEnvironmentPtr(), CPX_PARAM_BARALG, 3);
		CPXsetdblparam(cpx->getEnvironmentPtr(), CPX_PARAM_BAREPCOMP, 1.0e-6);
		if (par_->getIntParam("LOG_LEVEL") <= 5) {
			/** turn off all the messages (Barrier produce some even when loglevel = 0;) */
			si_->messageHandler()->setLogLevel(-1);
			CPXsetintparam(cpx->getEnvironmentPtr(), CPX_PARAM_BARDISPLAY, 0);
		} else {
			si_->messageHandler()->setLogLevel(1);
			CPXsetintparam(cpx->getEnvironmentPtr(), CPX_PARAM_BARDISPLAY, 1);
		}
	} else {
		si_->setHintParam(OsiDoPresolveInResolve, false);
		si_->setHintParam(OsiDoDualInResolve, false);
	}
#endif
	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

DSP_RTN_CODE DwMaster::solve() {
	BGN_TRY_CATCH

	itercnt_ = 0;
	t_total_ = CoinGetTimeOfDay();
	t_master_ = 0.0;
	t_colgen_ = 0.0;

	if (cols_generated_.size() == 0) {
		/** generate initial columns */
		DSP_RTN_CHECK_RTN_CODE(initialColumns());

		DSP_RTN_CHECK_RTN_CODE(solvePhase1());
	}

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
			if (primobj_ > feastol_)
				status_ = DSP_STAT_PRIM_INFEASIBLE;
			else {
				DSPdebugMessage("Converting to Phase 2.\n");
				DSP_RTN_CHECK_RTN_CODE(solvePhase2());
			}
		}
	}

	/** collect solutions */
	if (status_ == DSP_STAT_OPTIMAL ||
		status_ == DSP_STAT_FEASIBLE ||
		status_ == DSP_STAT_LIM_ITERorTIME) {

		/** recover original solution */
		CoinZeroN(primsol_, ncols_orig_);
		int cpos = 0;
		for (auto it = cols_generated_.begin(); it != cols_generated_.end(); it++) {
			if ((*it)->active_) {
				for (int i = 0; i < (*it)->x_.getNumElements(); ++i) {
					if ((*it)->x_.getIndices()[i] < ncols_orig_)
						primsol_[(*it)->x_.getIndices()[i]] += (*it)->x_.getElements()[i] * si_->getColSolution()[cpos];
				}
				cpos++;
			}
		}
#ifdef DSP_DEBUG_MORE
		printf("cols_generated_.size() %u, cpos %d,  ncols %d\n", cols_generated_.size(), cpos, si_->getNumCols());
		printf("lambda:\n");
		DspMessage::printArray(si_->getNumCols(), si_->getColSolution());
		//printf("x:\n");
		//DspMessage::printArray(ncols_orig_, primsol_);
		printf("Branch rows:\n");
		for (int i = nrows_core_; i < nrows_; ++i) {
			printf("primal solution [%d] %+e\n", branch_row_to_col_[i], primsol_[branch_row_to_col_[i]]);
			CoinShallowPackedVector row = si_->getMatrixByRow()->getVector(i);
			for (int k = 0, j = 0; k < row.getNumElements(); ++k) {
				if (j > 0 && j % 5 == 0) printf("\n");
				printf("  [%6d] %+e", row.getIndices()[k], row.getElements()[k]);
				j++;
			}
			printf("\n");
		}
#endif

#ifdef DSP_DEBUG
		for (int j = 0; j < ncols_orig_; ++j)
			if (primsol_[j] < node_clbd_[j] - 1.0e-6 || primsol_[j] > node_cubd_[j] + 1.0e-6) {
				status_ = DSP_STAT_UNKNOWN;
				printf("Col %d: primsol_ %e, clbd %e, cubd %e\n", j, primsol_[j], node_clbd_[j], node_cubd_[j]);
				throw CoinError("Primal solution violates the column bounds.", "solve", "DwMaster");
			}
#endif

#ifdef DSP_DEBUG_MORE
		for (int i = nrows_core_; i < nrows_; ++i) {
			CoinShallowPackedVector row = si_->getMatrixByRow()->getVector(i);
			for (int j = 0; j < row.getNumElements(); ++j)
				printf("col %d x [%e] lambda [%e]\n", row.getIndices()[j], row.getElements()[j], si_->getColSolution()[row.getIndices()[j]]);
			message_->print(0, "RowActivity[%d] %e in [%e %e] (ncols_inds %d)\n",
					i, si_->getRowActivity()[i], si_->getRowLower()[i], si_->getRowUpper()[i], branch_row_to_col_[i]);
		}

		for (int j = 0; j < ncols_orig_; ++j)
			if (primsol_[j] < node_clbd_[j] - 1.0e-6 || primsol_[j] > node_cubd_[j] + 1.0e-6) {
				message_->print(0, "Error: primsol_[%d] %e not in [%e %e].\n", j, primsol_[j], node_clbd_[j], node_cubd_[j]);
				CoinError("Solution violates bound constraints", "solve", "DwMaster");
				si_->writeMps("debug");
				return DSP_RTN_ERR;
			}
#endif
		/** heuristics */
		DSP_RTN_CHECK_RTN_CODE(heuristics());
	} else if (status_ == DSP_STAT_DUAL_INFEASIBLE) {
		//DSPdebug(si_->writeMps("master"));
	}

	/** always start with phase 2 */
	if (phase_ == 1) {
		/** delete auxiliary columns */
		si_->deleteCols(auxcolindices_.size(), &auxcolindices_[0]);
		auxcolindices_.clear();
		DSPdebugMessage("Phase 2 has %d rows and %d columns.\n", si_->getNumRows(), si_->getNumCols());

		/** set objective function coefficients */
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
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DwMaster::initialColumns() {
#define FREE_MEMORY     \
	FREE_ARRAY_PTR(piA)

	/** column generation info */
	double* piA = NULL;
	int nColsAdded;

	BGN_TRY_CATCH

	/** should consider phase 2 */
	phase_ = 2;

	/** allocate memory */
	std::vector<double> price(nrows_orig_ + nrows_conv_);
	piA = new double [ncols_orig_];

	/** initial price to generate columns */
	CoinZeroN(&price[0], nrows_orig_);
	CoinFillN(&price[0] + nrows_orig_, nrows_conv_, COIN_DBL_MAX);

	/** generate columns */
	DSP_RTN_CHECK_RTN_CODE(generateCols(&price[0], piA, dualobj_, nColsAdded));
	message_->print(0, "Generated %u initial columns. Initial dual bound %e\n", cols_generated_.size(), dualobj_);

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

		/** set zeros for the objective function coefficients */
		for (int j = 0; j < ncols; ++j)
			si_->setObjCoeff(j, 0.0);

		/** add auxiliary columns */
		double auxcolval;
		for (int i = 0; i < nrows_; ++i) {
			if (i < nrows_orig_ || i >= nrows_core_) {
				switch (si_->getRowSense()[i]) {
				case 'G':
					auxcolval = 1.0;
					si_->addCol(1, &i, &auxcolval, 0.0, COIN_DBL_MAX, 1.0);
					break;
				case 'L':
					auxcolval = -1.0;
					si_->addCol(1, &i, &auxcolval, 0.0, COIN_DBL_MAX, 1.0);
					break;
				case 'E':
				case 'R':
					auxcolval = 1.0;
					si_->addCol(1, &i, &auxcolval, 0.0, COIN_DBL_MAX, 1.0);
					auxcolval = -1.0;
					si_->addCol(1, &i, &auxcolval, 0.0, COIN_DBL_MAX, 1.0);
					break;
				default:
					break;
				}
			}
		}

		/** store indicies for the auxiliary columns */
		for (unsigned j = ncols; j < si_->getNumCols(); ++j)
			auxcolindices_.push_back(j);
		phase_ = 1;
		message_->print(3, "Phase 1 has %d rows and %d columns.\n", si_->getNumRows(), si_->getNumCols());
	}

	/** set parameters */
	worker_->setGapTolerance(0.0);
	//worker_->setTimeLimit(1.0e+20);

	DSP_RTN_CHECK_RTN_CODE(gutsOfSolve());

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

DSP_RTN_CODE DwMaster::solvePhase2() {
	BGN_TRY_CATCH

	if (phase_ == 1) {
		/** delete auxiliary columns */
		si_->deleteCols(auxcolindices_.size(), &auxcolindices_[0]);
		auxcolindices_.clear();
		message_->print(3, "Phase 2 has %d rows and %d columns.\n", si_->getNumRows(), si_->getNumCols());

		/** set objective function coefficients */
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

	/** set parameters */
	worker_->setGapTolerance(par_->getDblParam("DW/GAPTOL"));
	worker_->setTimeLimit(par_->getDblParam("DW/SUB/TIME_LIM"));

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

	double curLb = -COIN_DBL_MAX;

	/** column generation info */
	int nColsAdded = 0;

	/** timing results */
	double stime;

	BGN_TRY_CATCH

	/** allocate memory */
	price = new double [si_->getNumRows()];
	piA = new double [ncols_orig_];

	std::vector<double> prevsol(si_->getNumCols(), 0.0);

	/** use dual simplex after branching */
	if (!useBarrier_) {
		si_->setHintParam(OsiDoDualInResolve, true);
		si_->messageHandler()->setLogLevel(0);
	}

	/** solver the master problem */
	DSP_RTN_CHECK_RTN_CODE(solveMaster());

	/** use primal simplex after column generation */
	if (!useBarrier_)
		si_->setHintParam(OsiDoDualInResolve, false);

	while (status_ == DSP_STAT_OPTIMAL) {

		/** get price */
		CoinCopyN(si_->getRowPrice(), si_->getNumRows(), price);
#ifdef DSP_DEBUG_MORE
		printf("price:\n");
		DspMessage::printArray(si_->getNumRows(), price);
#endif

		/** column management */
		DSP_RTN_CHECK_RTN_CODE(reduceCols());

		/** generate columns */
		stime = CoinGetTimeOfDay();
		DSP_RTN_CHECK_RTN_CODE(generateCols(price, piA, curLb, nColsAdded));
		t_colgen_ += CoinGetTimeOfDay() - stime;

		/** subproblem solution may declare infeasibility. */
		if (status_ == DSP_STAT_PRIM_INFEASIBLE)
			break;

		/** update master */
		DSP_RTN_CHECK_RTN_CODE(updateModel(price, curLb));

		/** solver the master problem */
		DSP_RTN_CHECK_RTN_CODE(solveMaster());

		/** print information and increment iteration */
		printIterInfo();
		itercnt_++;

		/** termination test */
		if (terminationTest(nColsAdded))
			break;

#ifdef DSP_DEBUG_CPX
		char fname[128];
		sprintf(fname, "master%d", itercnt_);
		si_->writeMps(fname);
#endif
	}

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

DSP_RTN_CODE DwMaster::solveMaster() {
	BGN_TRY_CATCH

	if (useBarrier_) {
		OoqpEps* ooqp = dynamic_cast<OoqpEps*>(si_);
		if (ooqp) {
			if (phase_ == 1)
				ooqp->setOoqpStatus(0.0, -COIN_DBL_MAX, COIN_DBL_MAX);
			else {
				double epsilon = relgap_ / 1.1;
				if (primobj_ > 1.0e+20) epsilon = 0.0001;
				if (epsilon > 1.) epsilon = 1.;
				ooqp->setOoqpStatus(epsilon, dualobj_, bestprimobj_);
			}
		}
	}

	/** resolve */
	DSPdebugMessage("solve master problem (nrows %d, ncols %d).\n", si_->getNumRows(), si_->getNumCols());
	double stime = CoinGetTimeOfDay();
	si_->resolve();
	t_master_ += CoinGetTimeOfDay() - stime;
	convertCoinToDspStatus(si_, status_);
#if 0
	if (status_ == DSP_STAT_ABORT && useBarrier_) {
		message_->print(1, "Barrier solver detected numerical issues. Changed to simplex solver.\n");

		si_->setHintParam(OsiDoPresolveInResolve, false);
		si_->setHintParam(OsiDoDualInResolve, true);
		useBarrier_ = false;

		/** resolve */
		stime = CoinGetTimeOfDay();
		si_->resolve();
		t_master_ += CoinGetTimeOfDay() - stime;
		convertCoinToDspStatus(si_, status_);
	}
#endif

	/** calculate primal objective value */
	primobj_ = si_->getObjValue();
	absgap_ = fabs(primobj_-dualobj_);
	relgap_ = absgap_/(1.0e-10+fabs(primobj_));
	DSPdebugMessage("master status %d, primobj_ %e, absgap_ %e, relgap_ %e\n", status_, primobj_, absgap_, relgap_);

	/** store previous solution and basis */
	if (!useBarrier_) {
//		ws = dynamic_cast<CoinWarmStartBasis*>(si_->getWarmStart());
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

DSP_RTN_CODE DwMaster::restoreCols() {

	std::vector<int> existingCols;
	std::vector<int> delCols;

	BGN_TRY_CATCH

	/** delete auxiliary columns */
	si_->deleteCols(auxcolindices_.size(), &auxcolindices_[0]);
	auxcolindices_.clear();

	/** mark the existing columns to match with basis; and mark columns to delete */
	for (unsigned k = 0, j = 0; k < cols_generated_.size(); ++k)
		if (cols_generated_[k]->active_) {
			existingCols.push_back(k);
			delCols.push_back(j++);
		}

	if (delCols.size() > 0) {
		/** delete columns */
		si_->deleteCols(delCols.size(), &delCols[0]);

		/** add columns */
		for (unsigned k = 0; k < cols_generated_.size(); ++k) {
			cols_generated_[k]->active_ = true;
			si_->addCol(cols_generated_[k]->col_, cols_generated_[k]->lb_, cols_generated_[k]->ub_, cols_generated_[k]->obj_);
		}
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DwMaster::reduceCols() {
	BGN_TRY_CATCH
	/** FIXME:
	 * Be careful! Phase 2 can be infeasible by deleting columns.
	 * Do I really want to have this situation? */
	if (phase_ == 2) {
		std::vector<int> delcols;
		for (unsigned k = 0, j = 0; k < cols_generated_.size(); ++k) {
			if (cols_generated_[k]->active_) {
				/** age? */
				if (si_->getReducedCost()[j] < 1.0e-8)
					cols_generated_[k]->age_ = 0;
				else
					cols_generated_[k]->age_++;
				/** reduced cost fixing */
				if (dualobj_ + si_->getReducedCost()[j] - bestprimobj_ > -1.0e-10) {
					cols_generated_[k]->active_ = false;
					delcols.push_back(j);
				}
#if 0
				/** delete old cuts */
				if (cols_generated_[k]->active_ &&
						cols_generated_[k]->age_ >= par_->getIntParam("DW/MASTER/COL_AGE_LIM")) {
					cols_generated_[k]->active_ = false;
					delcols.push_back(j);
				}
#endif
				j++;
			}
		}
		if (delcols.size() > 0) {
			CoinWarmStartBasis* ws = NULL;
			if (useBarrier_ == false) {
				ws = dynamic_cast<CoinWarmStartBasis*>(si_->getWarmStart());
				ws->deleteColumns(delcols.size(), &delcols[0]);
			}
			si_->deleteCols(delcols.size(), &delcols[0]);
			if (useBarrier_ == false) {
				si_->setWarmStart(ws);
				FREE_PTR(ws)
			}
			message_->print(1, "Reduced cost fixing: removed %u columns.\n", delcols.size());

			si_->resolve();
			if (useBarrier_ == false)
				ws = dynamic_cast<CoinWarmStartBasis*>(si_->getWarmStart());
		}
	}
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
		double cutoff = price[nrows_orig_ + sind];
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

			/** convex combination constraints */
			if (statuses[s] != DSP_STAT_DUAL_INFEASIBLE)
				colvec.insert(nrows_orig_ + sind, 1.0);

			/** branching constraints */
			for (int i = 0; i < nrows_branch_; ++i) {
				int j = branch_row_to_col_[nrows_core_ + i];
				int sparse_index = x->findIndex(j);
				if (sparse_index == -1) continue;
				double val = x->getElements()[sparse_index];
				if (fabs(val) > 1.0e-10)
					colvec.insert(nrows_core_ + i, val);
			}

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
	DSPdebugMessage("Number of columns in the pool: %u\n", cols_generated_.size());

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
	lb = 0.0;
	double absprice = 0.0;
	const double* rlbd = si_->getRowLower();
	const double* rubd = si_->getRowUpper();
	for (int j = 0; j < nrows_; ++j) {
		if (j >= nrows_orig_ && j < nrows_core_)
			continue;
		absprice += fabs(price[j]);
		if (rlbd[j] > -1.0e+20)
			lb += price[j] * rlbd[j];
		else if (rubd[j] < 1.0e+20)
			lb += price[j] * rubd[j];
	}
	DSPdebugMessage("|price| %e, pia %e\n", absprice, lb);
	for (unsigned int s = 0; s < objs.size(); ++s) {
		lb += objs[s];
		DSPdebugMessage("subobj[%d] %e\n", s, objs[s]);
	}
	DSPdebugMessage("lb %e\n", lb);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

bool DwMaster::terminationTest(int nnewcols) {

	if (nnewcols == 0 || (!useBarrier_ && si_->getIterationCount() == 0))
		return true;

	bool term = false;

	if (phase_ == 1) {
		if (si_->isProvenOptimal() && primobj_ < feastol_) {
			status_ = DSP_STAT_FEASIBLE;
			DSPdebugMessage("Phase 1 found a feasible solution! %e < %e\n", si_->getObjValue(), feastol_);
			term = true;
		}
	}
	if (phase_ == 2) {
		if (iterlim_ <= itercnt_ ||
			relgap_ < par_->getDblParam("DW/GAPTOL") ||
			dualobj_ >= bestprimobj_) {
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
	for (int i = nrows_core_; i < nrows_; ++i)
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
		if (statuses[s] != DSP_STAT_OPTIMAL &&
				statuses[s] != DSP_STAT_DUAL_INFEASIBLE &&
				statuses[s] != DSP_STAT_LIM_ITERorTIME) {
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
	CoinCopyN(si_->getRowLower() + nrows_core_, nrows_branch_, rlbd);
	CoinCopyN(si_->getRowUpper() + nrows_core_, nrows_branch_, rubd);

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
		si_->setRowBounds(nrows_core_ + i, rlbd[i], rubd[i]);
		branchIndices.push_back(branch_row_to_col_[nrows_core_ + i]);
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
		if (org_obj_[branch_row_to_col_[nrows_core_ + i]] >= 0) continue;
		/** skip fixed bounds */
		if (rlbd[i] == rubd[i]) continue;
		si_->setRowBounds(nrows_core_ + i, 0.0, 0.0);
		branchIndices.push_back(branch_row_to_col_[nrows_core_ + i]);
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
			double x = si_->getRowActivity()[nrows_core_ + i];
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
					for (int i = 0; i < xlam.getNumElements(); ++i) {
						if (xlam.getIndices()[i] < ncols_orig_)
							bestprimsol_[xlam.getIndices()[i]] += xlam.getElements()[i];
					}
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
		int colind = branch_row_to_col_[nrows_core_ + i];
		if (org_obj_[colind] * direction <= 0) continue;

		double rounded = round(primsol_[colind]);
		si_->setRowBounds(nrows_core_ + i, rounded, rounded);
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
			double x = si_->getRowActivity()[nrows_core_ + i];
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
					for (int i = 0; i < xlam.getNumElements(); ++i) {
						if (xlam.getIndices()[i] < ncols_orig_)
							bestprimsol_[xlam.getIndices()[i]] += xlam.getElements()[i];
					}
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
				for (int i = 0; i < xlam.getNumElements(); ++i) {
					if (xlam.getIndices()[i] < ncols_orig_)
						primsol_[xlam.getIndices()[i]] += xlam.getElements()[i];
				}
				j++;
			}
		}

		/** find a fractional */
		findPhase = 0;
		branchIndex = -1;
		mindist = 1.0;
		isInteger = true;
		while (findPhase < 2 && branchIndex < 0) {
			for (int i = 0; i < nrows_branch_; ++i) {
				int colind = branch_row_to_col_[nrows_core_ + i];
				/** do not consider those with non-negative objective coefficients */
				if (findPhase == 0 && org_obj_[colind] >= 0) continue;
				/** do not consider those with negative objective coefficients */
				if (findPhase == 1 && org_obj_[colind] < 0) continue;

				/** skip if already fixed */
				if (si_->getRowUpper()[nrows_core_ + i] - si_->getRowLower()[nrows_core_ + i] < 1.0e-8)
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
			findPhase++;
		}

		/** found a fractional variable */
		if (branchIndex > -1) {
			/** keep the current node bounds */
			double rlbd_node = si_->getRowLower()[nrows_core_ + branchIndex];
			double rubd_node = si_->getRowUpper()[nrows_core_ + branchIndex];

			/** fix bounds */
			double rounded = round(branchValue);
			si_->setRowBounds(nrows_core_ + branchIndex, rounded, rounded);
			worker_->setColBounds(1, &branch_row_to_col_[nrows_core_ + branchIndex], &rounded, &rounded);
			message_->print(2, "Diving fixed variable %d [%e] to %e (discrepancy %u, depth %d).\n", branchIndex, branchValue, rounded, branchList.size(), depth);

			/** recursive call */
			status = status_;
			DSP_RTN_CHECK_RTN_CODE(gutsOfDive(branchList, depth+1));
			status_ = status;

			/** restore node bounds */
			si_->setRowBounds(nrows_core_ + branchIndex, rlbd_node, rubd_node);
			worker_->setColBounds(1, &branch_row_to_col_[nrows_core_ + branchIndex], &rlbd_node, &rubd_node);

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
			std::vector<int> branchIndices;
			std::vector<double> branchBounds;
			double* rlbd_tmp = NULL;
			double* rubd_tmp = NULL;
			if (nrows_branch_ > 0) {
				/** backup and fix all bounds */
				rlbd_tmp = new double [nrows_branch_];
				rubd_tmp = new double [nrows_branch_];
				CoinCopyN(si_->getRowLower() + nrows_core_, nrows_branch_, rlbd_tmp);
				CoinCopyN(si_->getRowUpper() + nrows_core_, nrows_branch_, rubd_tmp);
				for (int j = 0; j < nrows_branch_; ++j) {
					double rounded = round(primsol_[branch_row_to_col_[nrows_core_ + j]]);
					si_->setRowBounds(nrows_core_ + j, rounded, rounded);
					branchIndices.push_back(branch_row_to_col_[nrows_core_ + j]);
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
			}
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
					for (int i = 0; i < xlam.getNumElements(); ++i) {
						if (xlam.getIndices()[i] < ncols_orig_)
							primsol_[xlam.getIndices()[i]] += xlam.getElements()[i];
					}
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
			if (nrows_branch_ > 0) {
				for (int j = 0; j < nrows_branch_; ++j)
					si_->setRowBounds(nrows_core_ + j, rlbd_tmp[j], rubd_tmp[j]);
				worker_->setColBounds(branchIndices.size(), &branchIndices[0], rlbd_tmp, rubd_tmp);
				FREE_ARRAY_PTR(rlbd_tmp)
				FREE_ARRAY_PTR(rubd_tmp)
			}
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
		int row = nrows_core_ + i;
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
			for (int i = 0; i < xlam.getNumElements(); ++i) {
				if (xlam.getIndices()[i] < ncols_orig_)
					bestprimsol_[xlam.getIndices()[i]] += xlam.getElements()[i];
			}
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

	int findPhase = 0;
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
		branchingUp->bestBound_ = dualobj_;
		branchingDn->bestBound_ = dualobj_;
	} else {
		DSPdebugMessage("No branch object is found.\n");
	}

	END_TRY_CATCH_RTN(;,false)

	return branched;
}

void DwMaster::setBranchingObjects(const DspBranch* branchobj) {
	/** shouldn't be null */
	if (branchobj == NULL)
		return;

	std::vector<int> delrows;
	std::vector<int> delcols;

	/** adding columns */
	std::vector<int> col_inds;
	std::vector<double> col_elems;

	BGN_TRY_CATCH

	/** remove all the branching rows */
	if (nrows_branch_ > 0) {
		delrows.reserve(nrows_branch_);
		for (int i = nrows_core_; i < nrows_; ++i)
			delrows.push_back(i);
		si_->deleteRows(nrows_branch_, &delrows[0]);
		DSPdebugMessage("Deleted %d rows in the master.\n", nrows_branch_);
		nrows_branch_ = 0;
	}

	/** remove all columns */
	int ndelcols = si_->getNumCols();
	delcols.reserve(ndelcols);
	for (int j = 0; j < ndelcols; ++j)
		delcols.push_back(j);
	si_->deleteCols(ndelcols, &delcols[0]);
	DSPdebugMessage("Deleted %d columns in the master.\n", ndelcols);

	/** count nrows_branch_ */
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

	/** update number of rows */
	nrows_ = nrows_core_ + nrows_branch_;
	DSPdebugMessage("nrows_ %d nrows_core_ %d nrows_branch_ %d\n", nrows_, nrows_core_, nrows_branch_);

	/** reserve vector sizes */
	col_inds.reserve(nrows_);
	col_elems.reserve(nrows_);

	/** add branching rows */
	for (auto it = cols_generated_.begin(); it != cols_generated_.end(); it++) {
		if ((*it)->active_) {
			/** create a column for core rows */
			col_inds.clear();
			col_elems.clear();
			for (int i = 0; i < (*it)->col_.getNumElements(); ++i) {
				if ((*it)->col_.getIndices()[i] < nrows_core_) {
					col_inds.push_back((*it)->col_.getIndices()[i]);
					col_elems.push_back((*it)->col_.getElements()[i]);
				}
			}

			/** append column elements for the branching rows */
			for (unsigned j = 0, i = 0; j < branchobj->index_.size(); ++j) {
				int sparse_index = (*it)->x_.findIndex(branchobj->index_[j]);
				double val = 0.0;
				if (sparse_index > -1)
					val = (*it)->x_.getElements()[sparse_index];
				if (branchobj->lb_[j] > org_clbd_[branchobj->index_[j]]) {
					if (fabs(val) > 1.0e-10) {
						col_inds.push_back(nrows_core_+i);
						col_elems.push_back(val);
					}
					i++;
				}
				if (branchobj->ub_[j] < org_cubd_[branchobj->index_[j]]) {
					if (fabs(val) > 1.0e-10) {
						col_inds.push_back(nrows_core_+i);
						col_elems.push_back(val);
					}
					i++;
				}
			}

			/** assign the core-row column */
			(*it)->col_.clear();
			(*it)->col_.setVector(col_inds.size(), &col_inds[0], &col_elems[0]);

			/** add column */
			si_->addCol((*it)->col_, 0.0, COIN_DBL_MAX, (*it)->obj_);
		}
	}
	DSPdebugMessage("Appended dynamic columns in the master (%d / %u cols).\n", si_->getNumCols(), cols_generated_.size());

	/** restore column bounds */
	CoinCopyN(org_clbd_, ncols_orig_, node_clbd_);
	CoinCopyN(org_cubd_, ncols_orig_, node_cubd_);

	/** update column bounds at the current node */
	for (unsigned j = 0, irow = nrows_core_; j < branchobj->index_.size(); ++j) {
		node_clbd_[branchobj->index_[j]] = branchobj->lb_[j];
		node_cubd_[branchobj->index_[j]] = branchobj->ub_[j];
#ifdef DSP_DEBUG
		printf("Branch Obj: index %d lb %e ub %e\n", branchobj->index_[j], branchobj->lb_[j], branchobj->ub_[j]);
		CoinShallowPackedVector row = si_->getMatrixByRow()->getVector(irow);
		for (int i = 0, j = 0; i < row.getNumElements(); ++i) {
			if (j > 0 && j % 5 == 0) printf("\n");
			printf("  [%6d] %+e", row.getIndices()[i], row.getElements()[i]);
			j++;
		}
		printf("\n");
		irow++;
#endif
	}

	/** apply column bounds */
	std::vector<int> ncols_inds(ncols_orig_);
	CoinIotaN(&ncols_inds[0], ncols_orig_, 0);
	worker_->setColBounds(ncols_orig_, &ncols_inds[0], node_clbd_, node_cubd_);

	/** set known best bound */
	dualobj_ = branchobj->bestBound_;

	END_TRY_CATCH(;)
}

void DwMaster::printIterInfo() {
	message_->print(2, "[Phase %d] Iteration %3d: Master objective %e, ", phase_, itercnt_, primobj_);
	if (phase_ == 2) {
		if (dualobj_ > -1.0e+50)
			message_->print(2, "Lb %e (gap %.2f %%), ", dualobj_, relgap_*100);
		else
			message_->print(2, "Lb -Inf, ");
	}
	message_->print(2, "nrows %d, ncols %d, ", si_->getNumRows(), si_->getNumCols());
	if (!useBarrier_)
		message_->print(2, "itercnt %d, ", si_->getIterationCount());
	message_->print(2, "timing (total %.2f, master %.2f, gencols %.2f), statue %d\n",
			CoinGetTimeOfDay() - t_total_, t_master_, t_colgen_, status_);
}
