//
// Created by Kibaek Kim on 8/27/16.
//

// #define DSP_DEBUG

#include "SolverInterface/DspOsiClp.h"
#include "SolverInterface/DspOsiCpx.h"
#include "SolverInterface/DspOsiGrb.h"
#include "Model/TssModel.h"
//#include "SolverInterface/OoqpEps.h"
#include "Solver/DantzigWolfe/DwMaster.h"
#include "Utility/DspUtility.h"

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
mat_orig_(NULL),
branchObj_(NULL),
itercnt_(0),
ngenerated_(0),
t_start_(0.0),
t_total_(0.0),
t_master_(0.0),
t_colgen_(0.0) {
	useBarrier_ = par_->getBoolParam("DW/MASTER/IPM");
}

DwMaster::DwMaster(const DwMaster& rhs):
DecSolver(rhs),
useBarrier_(rhs.useBarrier_),
phase_(rhs.phase_),
auxcolindices_(rhs.auxcolindices_),
worker_(rhs.worker_),
ncols_orig_(rhs.ncols_orig_),
ncols_start_(rhs.ncols_start_),
nrows_(rhs.nrows_),
nrows_orig_(rhs.nrows_orig_),
nrows_conv_(rhs.nrows_conv_),
nrows_core_(rhs.nrows_core_),
nrows_branch_(rhs.nrows_branch_),
#ifdef USE_ROW_TO_COL
branch_row_to_col_(rhs.branch_row_to_col_),
#else
branch_row_to_vec_(rhs.branch_row_to_vec_),
#endif
clbd_orig_(rhs.clbd_orig_),
cubd_orig_(rhs.cubd_orig_),
obj_orig_(rhs.obj_orig_),
ctype_orig_(rhs.ctype_orig_),
rlbd_orig_(rhs.rlbd_orig_),
rubd_orig_(rhs.rubd_orig_),
clbd_node_(rhs.clbd_node_),
cubd_node_(rhs.cubd_node_),
branchObj_(rhs.branchObj_),
bestprimsol_orig_(rhs.bestprimsol_orig_),
itercnt_(rhs.itercnt_),
ngenerated_(rhs.ngenerated_),
log_time_(rhs.log_time_),
log_bestdual_bounds_(rhs.log_bestdual_bounds_),
log_bestprim_bounds_(rhs.log_bestprim_bounds_),
t_start_(rhs.t_start_),
t_total_(rhs.t_total_),
t_master_(rhs.t_master_),
t_colgen_(rhs.t_colgen_),
status_subs_(rhs.status_subs_) {
	mat_orig_ = new CoinPackedMatrix(*(rhs.mat_orig_));
	for (auto it = rhs.cols_generated_.begin(); it != rhs.cols_generated_.end(); it++)
		cols_generated_.push_back(new DwCol(**it));
	for (auto it = rhs.recent_subsols_.begin(); it != rhs.recent_subsols_.end(); it++)
		recent_subsols_.push_back(new CoinPackedVector(**it));
	for (auto it = rhs.stored_solutions_.begin(); it != rhs.stored_solutions_.end(); it++)
		stored_solutions_.push_back(new CoinPackedVector(**it));
}

/** copy operator */
DwMaster& DwMaster::operator=(const DwMaster& rhs) {
	DecSolver::operator=(rhs);
	osi_ = rhs.osi_;
	useBarrier_ = rhs.useBarrier_;
	phase_ = rhs.phase_;
	auxcolindices_ = rhs.auxcolindices_;
	worker_ = rhs.worker_;
	ncols_orig_ = rhs.ncols_orig_;
	ncols_start_ = rhs.ncols_start_;
	nrows_ = rhs.nrows_;
	nrows_orig_ = rhs.nrows_orig_;
	nrows_conv_ = rhs.nrows_conv_;
	nrows_core_ = rhs.nrows_core_;
	nrows_branch_ = rhs.nrows_branch_;
#ifdef USE_ROW_TO_COL
	branch_row_to_col_ = rhs.branch_row_to_col_;
#else
	branch_row_to_vec_ = rhs.branch_row_to_vec_;
#endif
	clbd_orig_ = rhs.clbd_orig_;
	cubd_orig_ = rhs.cubd_orig_;
	obj_orig_ = rhs.obj_orig_;
	ctype_orig_ = rhs.ctype_orig_;
	rlbd_orig_ = rhs.rlbd_orig_;
	rubd_orig_ = rhs.rubd_orig_;
	clbd_node_ = rhs.clbd_node_;
	cubd_node_ = rhs.cubd_node_;
	branchObj_ = rhs.branchObj_;
	itercnt_ = rhs.itercnt_;
	ngenerated_ = rhs.ngenerated_;
	t_start_ = rhs.t_start_;
	t_total_ = rhs.t_total_;
	t_master_ = rhs.t_master_;
	t_colgen_ = rhs.t_colgen_;
	status_subs_ = rhs.status_subs_;
	mat_orig_ = new CoinPackedMatrix(*(rhs.mat_orig_));
	for (auto it = rhs.cols_generated_.begin(); it != rhs.cols_generated_.end(); it++)
		cols_generated_.push_back(new DwCol(**it));
	for (auto it = rhs.stored_solutions_.begin(); it != rhs.stored_solutions_.end(); it++)
		stored_solutions_.push_back(new CoinPackedVector(**it));
	for (auto it = rhs.recent_subsols_.begin(); it != rhs.recent_subsols_.end(); it++)
		recent_subsols_.push_back(new CoinPackedVector(**it));
	return *this;
}

DwMaster::~DwMaster() {
	osi_ = NULL;
	branchObj_ = NULL;
	FREE_PTR(mat_orig_);
	for (unsigned i = 0; i < cols_generated_.size(); ++i)
		FREE_PTR(cols_generated_[i]);
	for (unsigned i = 0; i < stored_solutions_.size(); ++i)
		FREE_PTR(stored_solutions_[i]);
	for (unsigned i = 0; i < recent_subsols_.size(); ++i)
		FREE_PTR(recent_subsols_[i]);
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
		DSPdebug(org_mat->verifyMtx(4));

		int nscen = tss->getNumScenarios();
		int ncols_first_stage = tss->getNumCols(0);
		int ncols = org_mat->getNumCols() + ncols_first_stage * (nscen - 1);
		const double* probability = tss->getProbability();
		DSPdebugMessage("nscen %d ncols_first_stage %d ncols %d\n", nscen, ncols_first_stage, ncols);

		mat_orig_ = new CoinPackedMatrix(org_mat->isColOrdered(), 0, 0);
		mat_orig_->setDimensions(0, ncols);

		/** add non-anticipativity constraints in the following form:
		   x_1 - x_2      = 0
		        x_2 - x_3 = 0 
		  -x_1      + x_3 = 0
		*/
		int indices[2];
		double elements[] = {1.0, -1.0};
		for (int i = 0; i < nscen; ++i) {
			if (i < nscen - 1) {
				for (int j = 0; j < ncols_first_stage; ++j) {
					indices[0] = i * ncols_first_stage + j;
					indices[1] = (i+1) * ncols_first_stage + j;
					mat_orig_->appendRow(2, indices, elements);
				}
			} else {
				for (int j = 0; j < ncols_first_stage; ++j) {
					indices[0] = i * ncols_first_stage + j;
					indices[1] = j;
					mat_orig_->appendRow(2, indices, elements);
				}
			}
		}
		DSPdebug(mat_orig_->verifyMtx(4));

		clbd_orig_.resize(ncols);
		cubd_orig_.resize(ncols);
		ctype_orig_.resize(ncols);
		obj_orig_.resize(ncols);
		rlbd_orig_.resize(mat_orig_->getNumRows());
		rubd_orig_.resize(mat_orig_->getNumRows());
		for (int s = 0; s < nscen; ++s) {
			std::copy(org_ctype, org_ctype + ncols_first_stage, ctype_orig_.begin() + s * ncols_first_stage);
			std::copy(org_clbd, org_clbd + ncols_first_stage, clbd_orig_.begin() + s * ncols_first_stage);
			std::copy(org_cubd, org_cubd + ncols_first_stage, cubd_orig_.begin() + s * ncols_first_stage);
			for (int j = 0; j < ncols_first_stage; ++j)
				obj_orig_[s * ncols_first_stage + j] = org_obj[j] * probability[s];
		}
		std::copy(org_ctype + ncols_first_stage, org_ctype + ncols - (nscen-1) * ncols_first_stage, ctype_orig_.begin() + nscen * ncols_first_stage);
		std::copy(org_clbd + ncols_first_stage, org_clbd + ncols - (nscen-1) * ncols_first_stage, clbd_orig_.begin() + nscen * ncols_first_stage);
		std::copy(org_cubd + ncols_first_stage, org_cubd + ncols - (nscen-1) * ncols_first_stage, cubd_orig_.begin() + nscen * ncols_first_stage);
		std::fill(obj_orig_.begin() + nscen * ncols_first_stage, obj_orig_.end(), 0.0);
		std::fill(rlbd_orig_.begin(), rlbd_orig_.end(), 0.0);
		std::fill(rubd_orig_.begin(), rubd_orig_.end(), 0.0);
	} else {
		/** retrieve the original master problem structure */
		model_->decompose(0, NULL, 0, NULL, NULL, NULL,
				mat_orig_, org_clbd, org_cubd, org_ctype, org_obj, org_rlbd, org_rubd);
		clbd_orig_.assign(org_clbd, org_clbd + mat_orig_->getNumCols());
		cubd_orig_.assign(org_cubd, org_cubd + mat_orig_->getNumCols());
		ctype_orig_.assign(org_ctype, org_ctype + mat_orig_->getNumCols());
		obj_orig_.assign(org_obj, org_obj + mat_orig_->getNumCols());
		rlbd_orig_.assign(org_rlbd, org_rlbd + mat_orig_->getNumRows());
		rubd_orig_.assign(org_rubd, org_rubd + mat_orig_->getNumRows());
	}

	ncols_orig_ = mat_orig_->getNumCols(); /**< number of columns in the original master */
	nrows_orig_ = mat_orig_->getNumRows(); /**< number of rows in the original master */
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
	BGN_TRY_CATCH

	/** allocate memory */
	std::vector<double> rlbd(nrows_), rubd(nrows_);
	clbd_node_ = clbd_orig_;
	cubd_node_ = cubd_orig_;

	/** create column-wise matrix and set number of rows */
	std::shared_ptr<CoinPackedMatrix> mat(new CoinPackedMatrix(true, 0, 0));
	mat->setDimensions(nrows_, 0);

	/** Set row bounds */
	std::fill(rlbd.begin(), rlbd.begin() + nrows_conv_, 1.0);
	std::fill(rubd.begin(), rubd.begin() + nrows_conv_, 1.0);
	std::copy(rlbd_orig_.begin(), rlbd_orig_.end(), rlbd.begin() + nrows_conv_);
	std::copy(rubd_orig_.begin(), rubd_orig_.end(), rubd.begin() + nrows_conv_);

	/** create solver */
	osi_ = createDspOsi();
	if (!osi_) throw CoinError("Failed to create DspOsi", "createProblem", "DwMaster");

	/** set log level */
	osi_->setLogLevel(par_->getIntParam("DW/MASTER/SOLVER/LOG_LEVEL"));

	/** load problem data */
	getSiPtr()->loadProblem(*(mat.get()), NULL, NULL, NULL, &rlbd[0], &rubd[0]);

	dualsol_.reserve(getSiPtr()->getNumRows());

	DSP_RTN_CHECK_RTN_CODE(initialOsiSolver());

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DwMaster::initialOsiSolver() {
	BGN_TRY_CATCH
	/** set hint parameters */
	if (useBarrier_) {
		//
	} else {
		osi_->setNumCores(1);
		getSiPtr()->setHintParam(OsiDoPresolveInResolve, false);
		getSiPtr()->setHintParam(OsiDoDualInResolve, false);
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
			osi_->setLogLevel(-1);
			CPXsetintparam(cpx->getEnvironmentPtr(), CPX_PARAM_BARDISPLAY, 0);
		} else {
			osi_->setLogLevel(1);
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
	t_start_ = CoinGetTimeOfDay();
	t_total_ = 0.0;
	t_master_ = 0.0;
	t_colgen_ = 0.0;

	if (cols_generated_.size() == 0) {
		/** generate initial columns */
		DSP_RTN_CHECK_RTN_CODE(initialColumns());
		bestdualobj_ = std::max(bestdualobj_, dualobj_);

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

	/** switch to phase 2 */
	DSP_RTN_CHECK_RTN_CODE(switchToPhase2());

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DwMaster::initialColumns() {
	BGN_TRY_CATCH

	/** should consider phase 2 */
	phase_ = 2;

	/** allocate memory */
	dualsol_.resize(nrows_conv_ + nrows_orig_);

	/** initial price to generate columns */
	CoinFillN(&dualsol_[0], nrows_conv_, COIN_DBL_MAX);
	CoinZeroN(&dualsol_[nrows_conv_], nrows_orig_);

	/** generate columns */
	DSP_RTN_CHECK_RTN_CODE(generateCols());
	message_->print(3, "Generated %u initial columns. Initial dual bound %e\n", ngenerated_, dualobj_);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

DSP_RTN_CODE DwMaster::solvePhase1() {
	BGN_TRY_CATCH

	if (phase_ == 2) {
		/** set phase 1 problem */
		int ncols = getSiPtr()->getNumCols();

		/** set zeros for the objective function coefficients */
		for (int j = 0; j < ncols; ++j)
			getSiPtr()->setObjCoeff(j, 0.0);

		/** add auxiliary columns */
		double auxcolval;
		for (int i = nrows_conv_; i < nrows_; ++i) {
			switch (getSiPtr()->getRowSense()[i]) {
			case 'G':
				auxcolval = 1.0;
				getSiPtr()->addCol(1, &i, &auxcolval, 0.0, COIN_DBL_MAX, 1.0);
				break;
			case 'L':
				auxcolval = -1.0;
				getSiPtr()->addCol(1, &i, &auxcolval, 0.0, COIN_DBL_MAX, 1.0);
				break;
			case 'E':
			case 'R':
				auxcolval = 1.0;
				getSiPtr()->addCol(1, &i, &auxcolval, 0.0, COIN_DBL_MAX, 1.0);
				auxcolval = -1.0;
				getSiPtr()->addCol(1, &i, &auxcolval, 0.0, COIN_DBL_MAX, 1.0);
				break;
			default:
				break;
			}
		}

		/** store indicies for the auxiliary columns */
		for (unsigned j = ncols; j < getSiPtr()->getNumCols(); ++j)
			auxcolindices_.push_back(j);
		phase_ = 1;
		message_->print(3, "Phase 1 has %d rows and %d columns.\n", getSiPtr()->getNumRows(), getSiPtr()->getNumCols());
	}

	/** set parameters */
	worker_->setGapTolerance(0.0);
	//worker_->setTimeLimit(1.0e+20);

	/** use dual simplex after branching */
	if (!useBarrier_) {
		getSiPtr()->setHintParam(OsiDoDualInResolve, true);
		osi_->setLogLevel(0);
	}

	DSP_RTN_CHECK_RTN_CODE(gutsOfSolve());

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

DSP_RTN_CODE DwMaster::solvePhase2() {
	BGN_TRY_CATCH

	/** switch to phase 2 */
	DSP_RTN_CHECK_RTN_CODE(switchToPhase2());

	/** set parameters */
	worker_->setGapTolerance(par_->getDblParam("DW/GAPTOL"));
	worker_->setTimeLimit(par_->getDblParam("DW/SUB/TIME_LIM"));

	/** use dual simplex after branching */
	if (!useBarrier_) {
		getSiPtr()->setHintParam(OsiDoDualInResolve, true);
		osi_->setLogLevel(0);
	}

	DSP_RTN_CHECK_RTN_CODE(gutsOfSolve());
	//DSPdebug(si_->writeMps("afterPhase2"));

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

DSP_RTN_CODE DwMaster::gutsOfSolve() {

	double stime; /**< timing results */
	int prev_nsols;

	BGN_TRY_CATCH

	dualobj_ = -COIN_DBL_MAX;

	/** solver the master problem */
	DSP_RTN_CHECK_RTN_CODE(solveMaster());

	/** print information and increment iteration */
	printIterInfo();
	itercnt_++;

	int num_removed = 0, num_restored = 0;
	while (status_ == DSP_STAT_OPTIMAL && terminationTest() == false) {

		/** column management */
		num_removed = 0;
		DSP_RTN_CHECK_RTN_CODE(reduceCols(num_removed));

		if (num_removed > 0)
			message_->print(1, "Removed %u columns.\n", num_removed);
		else {
			/** restore columns */
			num_restored = 0;
			DSP_RTN_CHECK_RTN_CODE(restoreCols(num_restored));

			if (num_restored > 0)
				message_->print(1, "Restored %u columns.\n", num_restored);
			else {
				stime = CoinGetTimeOfDay();
				/** generate columns */
				prev_nsols = stored_solutions_.size();
				DSP_RTN_CHECK_RTN_CODE(generateCols());

				/** generate extra columns by fixing the first-stage variables for SMIP */
				DSP_RTN_CHECK_RTN_CODE(generateColsByFix(stored_solutions_.size() - prev_nsols));

				t_colgen_ += CoinGetTimeOfDay() - stime;

				/** subproblem solution may declare infeasibility. */
				for (auto st = status_subs_.begin(); st != status_subs_.end(); st++) {
					//DSPdebugMessage("subproblem status %d\n", *st);
					if (*st == DSP_STAT_PRIM_INFEASIBLE) {
						status_ = DSP_STAT_PRIM_INFEASIBLE;
						break;
					}
				}
				if (status_ == DSP_STAT_PRIM_INFEASIBLE)
					break;

				// check time limit
				if (time_remains_ <= (CoinGetTimeOfDay() - t_start_)) {
					status_ = DSP_STAT_STOPPED_TIME;
					break;
				}

				/** update master */
				DSP_RTN_CHECK_RTN_CODE(updateModel());
			}
		}

		/** solver the master problem */
		stime = CoinGetTimeOfDay();
		DSP_RTN_CHECK_RTN_CODE(solveMaster());
		t_master_ += CoinGetTimeOfDay() - stime;
		t_total_ = CoinGetTimeOfDay() - t_start_;

		/** print information and increment iteration */
		printIterInfo();
		itercnt_++;

#ifdef DSP_DEBUG
		char fname[128];
		sprintf(fname, "master%d", itercnt_);
		osi_->si_->writeMps(fname);
#endif
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DwMaster::solveMaster() {
	BGN_TRY_CATCH

	if (useBarrier_) {
		/** TODO: Replace this by PIPS */
		/*
		OoqpEps* ooqp = dynamic_cast<OoqpEps*>(si_);
		if (ooqp) {
			if (phase_ == 1)
				ooqp->setOoqpStatus(0.0, -COIN_DBL_MAX, COIN_DBL_MAX);
			else {
				double epsilon = getRelApproxGap() / 1.1;
				if (primobj_ > 1.0e+20) epsilon = 0.0001;
				if (epsilon > 1.) epsilon = 1.;
				ooqp->setOoqpStatus(epsilon, bestdualobj_, bestprimobj_);
			}
		}
		*/
	} else {
		getSiPtr()->setHintParam(OsiDoDualInResolve, false);
	}

	/** resolve */
	DSPdebugMessage("solve master problem (nrows %d, ncols %d).\n", getSiPtr()->getNumRows(), getSiPtr()->getNumCols());
	getSiPtr()->resolve();
	status_ = osi_->status();
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
		convertOsiToDspStatus(si_, status_);
	}
#endif

	/** calculate primal objective value */
	primobj_ = osi_->getPrimObjValue();
	DSPdebugMessage("master status %d, primobj_ %e, absgap_ %e, relgap_ %e\n", status_, primobj_, getAbsApproxGap(), getRelApproxGap());

	/** retrieve price */
	dualsol_.assign(getSiPtr()->getRowPrice(), getSiPtr()->getRowPrice() + getSiPtr()->getNumRows());

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

DSP_RTN_CODE DwMaster::updateCol(DwCol* col) {
	BGN_TRY_CATCH

	std::vector<int> col_inds;
	std::vector<double> col_elems;
	col_inds.reserve(nrows_);
	col_elems.reserve(nrows_);

	/** create a column for core rows */
	for (int i = 0; i < col->col_.getNumElements(); ++i) {
		if (col->col_.getIndices()[i] < nrows_core_) {
			col_inds.push_back(col->col_.getIndices()[i]);
			col_elems.push_back(col->col_.getElements()[i]);
		}
	}

#ifdef USE_ROW_TO_COL
	for (int i = 0; i < nrows_branch_; ++i) {
		int j = branch_row_to_col_[nrows_core_ + i];
		int sparse_index = col->x_.findIndex(j);
		if (sparse_index == -1) continue;
		double val = col->x_.getElements()[sparse_index];
		if (fabs(val) > 1.0e-10) {
			col_inds.push_back(nrows_core_ + i);
			col_elems.push_back(val);
		}
	}
#else
	if (nrows_branch_ > 0) {
		double* dx = col->x_.denseVector(ncols_orig_);
		for (int i = 0; i < nrows_branch_; ++i) {
			double val = branch_row_to_vec_[nrows_core_+i].dotProduct(dx);
			if (fabs(val) > 1.0e-10) {
				col_inds.push_back(nrows_core_ + i);
				col_elems.push_back(val);
			}
		}
		FREE_ARRAY_PTR(dx);
	}
#endif

	/** assign the core-row column */
	col->col_.setVector(col_inds.size(), &col_inds[0], &col_elems[0]);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

DSP_RTN_CODE DwMaster::restoreCols(int &num_restored) {
	BGN_TRY_CATCH

	num_restored = 0;

	/** do only at phase 2 */
	if (phase_ == 2) {
		/** mark the existing columns to match with basis; and mark columns to delete */
		for (unsigned k = 0; k < cols_generated_.size(); ++k) {
			if (cols_generated_[k]->active_) continue;
			if (cols_generated_[k]->age_ < COIN_INT_MAX) {
				/** update column */
				DSP_RTN_CHECK_RTN_CODE(updateCol(cols_generated_[k]));

				/** add column */
				double lhs = cols_generated_[k]->col_.dotProduct(&dualsol_[0]);
				if (lhs - cols_generated_[k]->obj_ > 1.0e-6) {
					cols_generated_[k]->active_ = true;
					cols_generated_[k]->age_ = 0;

					/** set master index */
					cols_generated_[k]->master_index_ = getSiPtr()->getNumCols();

					/** add column */
					getSiPtr()->addCol(cols_generated_[k]->col_.getNumElements(), 
						cols_generated_[k]->col_.getIndices(),
						cols_generated_[k]->col_.getElements(),
						cols_generated_[k]->lb_, cols_generated_[k]->ub_, cols_generated_[k]->obj_);
					num_restored++;
				}
			}
		}
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

DSP_RTN_CODE DwMaster::reduceCols(int &num_removed) {
	BGN_TRY_CATCH
	/** FIXME:
	 * Be careful! Phase 2 can be infeasible by deleting columns.
	 * Do I really want to have this situation? */
	if (phase_ == 2) {
		std::vector<int> delcols;
		for (unsigned k = 0, j = 0; k < cols_generated_.size(); ++k) {
			if (cols_generated_[k]->active_) {
				/** age? */
				if (getSiPtr()->getReducedCost()[j] < 1.0e-8)
					cols_generated_[k]->age_ = 0;
				else
					cols_generated_[k]->age_++;
				/** reduced cost fixing */
				if (bestdualobj_ + getSiPtr()->getReducedCost()[j] - bestprimobj_ > -1.0e-10) {
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
				ws = dynamic_cast<CoinWarmStartBasis*>(getSiPtr()->getWarmStart());
				ws->deleteColumns(delcols.size(), &delcols[0]);
			}
			getSiPtr()->deleteCols(delcols.size(), &delcols[0]);
			if (useBarrier_ == false) {
				getSiPtr()->setWarmStart(ws);
				FREE_PTR(ws)
			}
			message_->print(1, "Reduced cost fixing: removed %u columns.\n", delcols.size());

			getSiPtr()->resolve();
			if (useBarrier_ == false)
				ws = dynamic_cast<CoinWarmStartBasis*>(getSiPtr()->getWarmStart());
			
		}

		num_removed = delcols.size();
	}
	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

DSP_RTN_CODE DwMaster::generateCols() {
	std::vector<double> piA;
	/** column generation info */
	std::vector<int> subinds;
	std::vector<double> subcxs;
	std::vector<double> subobjs;
	std::vector<CoinPackedVector*> subsols;

	BGN_TRY_CATCH

	/** calculate pi^T A */
	DSP_RTN_CHECK_RTN_CODE(calculatePiA(piA));

	// set time limit
	double sub_timlim = par_->getDblParam("DW/SUB/TIME_LIM");
	worker_->setTimeLimit(CoinMin(sub_timlim, time_remains_ - (CoinGetTimeOfDay() - t_start_)));

	/** generate columns */
	DSP_RTN_CHECK_RTN_CODE(
			worker_->generateCols(phase_, &piA[0], subinds, status_subs_, subcxs, subobjs, subsols));
	DSPdebugMessage("status_subs_.size() %lu\n", status_subs_.size());

	// reset time limit
	par_->setDblParam("DW/SUB/TIME_LIM", sub_timlim);

	/** any subproblem primal/dual infeasible? */
	bool isInfeasible = false;
	for (auto status = status_subs_.begin(); status != status_subs_.end(); status++)
		if (*status == DSP_STAT_PRIM_INFEASIBLE ||
			*status == DSP_STAT_DUAL_INFEASIBLE) {
			isInfeasible = true;
			break;
		}

	if (!isInfeasible) {
		/** calculate lower bound */
		if (phase_ == 2) {
			DSP_RTN_CHECK_RTN_CODE(getLagrangianBound(subobjs));
			DSPdebugMessage("Current lower bound %e, best lower bound %e\n", dualobj_, bestdualobj_);
		}

		/** clear recent subproblem solutions */
		for (unsigned i = 0; i < recent_subsols_.size(); ++i)
			FREE_PTR(recent_subsols_[i]);
		recent_subsols_.clear();

		/** assign the current subproblem solutions */
		recent_subsols_.reserve(subsols.size());
		for (unsigned i = 0; i < subsols.size(); ++i)
			recent_subsols_.push_back(subsols[i]);

		/** create and add columns */
		DSP_RTN_CHECK_RTN_CODE(
				addCols(subinds, status_subs_, subcxs, subobjs, subsols));

		if (model_->isStochastic() && 
			par_->getIntParam("DW/MAX_EVAL_UB") > 0 &&
			par_->getIntParam("DW/EVAL_UB") >= 0) {
			/** maximum number of solutions to evaluate */
			int max_stores = par_->getIntParam("DW/MAX_EVAL_UB");

			/** store solutions to distribute */
			TssModel* tss = dynamic_cast<TssModel*>(model_);
			for (unsigned i = 0; i < subsols.size(); ++i) {
				if (max_stores == 0) break;
				/** assign first-stage solution value */
				CoinPackedVector* first_stage_solution = new CoinPackedVector;
				first_stage_solution->reserve(tss->getNumCols(0));
				for (int j = 0; j < subsols[i]->getNumElements(); ++j) {
					if (subsols[i]->getIndices()[j] >= tss->getNumScenarios() * tss->getNumCols(0))
						break;
					first_stage_solution->insert(subsols[i]->getIndices()[j] % tss->getNumCols(0), subsols[i]->getElements()[j]);
					// printf("i %d first_stage_solution[%d] = %e\n", i, subsols[i]->getIndices()[j] % tss->getNumCols(0), subsols[i]->getElements()[j]);
				}
				/** store it if not duplicated */
				if (!duplicateVector(first_stage_solution, stored_solutions_)) {
					stored_solutions_.push_back(first_stage_solution);
					/** count */
					max_stores--;
				}
			}
		}

	} else
		ngenerated_ = 0;

	/** free memory for subproblem solutions */
	for (unsigned i = 0; i < subsols.size(); ++i)
		subsols[i] = NULL;
	subsols.clear();

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

DSP_RTN_CODE DwMaster::generateColsByFix(
		int nsols /**< [in] number of solutions to evaluate in LIFO way */) {

	if (model_->isStochastic() == false) 
		return DSP_RTN_OK;
	if (nsols <= 0)
		return DSP_RTN_OK;

	/** column generation info */
	std::vector<int> subinds;
	std::vector<int> substatuses;
	std::vector<double> subobjs;
	std::vector<CoinPackedVector*> subsols;
	std::vector<CoinPackedVector*> solutions_to_evaluate(nsols, NULL);

	BGN_TRY_CATCH

	TssModel* tss = dynamic_cast<TssModel*>(model_);

	// set time limit
	double sub_timlim = par_->getDblParam("DW/SUB/TIME_LIM");

	for (unsigned i = stored_solutions_.size() - 1, j = 0; i >= stored_solutions_.size() - nsols; --i)
		solutions_to_evaluate[j++] = stored_solutions_[i];

	std::vector<double> first_stage_solution_dense(tss->getNumCols(0));
	for (unsigned i = 0; i < solutions_to_evaluate.size(); ++i) {
		/** generate columns */
		std::fill(first_stage_solution_dense.begin(), first_stage_solution_dense.end(), 0.0);
		for (int j = 0; j < solutions_to_evaluate[i]->getNumElements(); ++j)
			first_stage_solution_dense[solutions_to_evaluate[i]->getIndices()[j]] = solutions_to_evaluate[i]->getElements()[j];

		worker_->setTimeLimit(CoinMin(sub_timlim, time_remains_ - (CoinGetTimeOfDay() - t_start_)));
		DSP_RTN_CHECK_RTN_CODE(
				worker_->generateColsByFix(&first_stage_solution_dense[0], subinds, substatuses, subobjs, subsols));

		/** any subproblem primal/dual infeasible? */
		bool isInfeasibleFix = false;
		for (auto status = substatuses.begin(); status != substatuses.end(); status++)
			if (*status == DSP_STAT_PRIM_INFEASIBLE ||
				*status == DSP_STAT_DUAL_INFEASIBLE) {
				isInfeasibleFix = true;
				break;
			}

		if (!isInfeasibleFix) {
			/** create and add columns */
			DSP_RTN_CHECK_RTN_CODE(
					addCols(subinds, substatuses, subobjs, subobjs, subsols));

			double newbound = 0.0;
			for (auto it = subobjs.begin(); it != subobjs.end(); it++)
				newbound += *it;
			if (newbound < bestprimobj_) {
				message_->print(3, "  Found new primal bound %e (< %e)\n", newbound, bestprimobj_);
				bestprimobj_ = newbound;
				bestprimsol_orig_.assign(ncols_orig_, 0.0);
				for (unsigned s = 0; s < subsols.size(); s++)
					for (int j = 0; j < subsols[s]->getNumElements(); ++j)
						bestprimsol_orig_[subsols[s]->getIndices()[j]] = subsols[s]->getElements()[j];
			}
		}

		/** clean up */
		subinds.clear();
		substatuses.clear();
		subobjs.clear();
		for (unsigned i = 0; i < subsols.size(); ++i)
			FREE_PTR(subsols[i]);
		subsols.clear();
	}

	// reset time limit
	par_->setDblParam("DW/SUB/TIME_LIM", sub_timlim);

	/** clear solution vector */
	for (unsigned i = 0; i < solutions_to_evaluate.size(); ++i)
		solutions_to_evaluate[i] = NULL;
	solutions_to_evaluate.clear();

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

DSP_RTN_CODE DwMaster::addCols(
		std::vector<int>& indices,           /**< [in] subproblem indices corresponding to cols*/
		std::vector<int>& statuses,          /**< [in] subproblem solution status */
		std::vector<double>& cxs,            /**< [in] solution times original objective coefficients */
		std::vector<double>& objs,           /**< [in] subproblem objective values */
		std::vector<CoinPackedVector*>& sols /**< [in] subproblem solutions */) {
#define FREE_MEMORY \
	FREE_ARRAY_PTR(Ax)

	double* Ax = NULL;
	CoinPackedVector colvec;

	BGN_TRY_CATCH

	/** allocate memory */
	Ax = new double [nrows_orig_];
	colvec.reserve(nrows_);

	/** reset counter */
	ngenerated_ = 0;

	for (unsigned int s = 0; s < indices.size(); ++s) {
		int sind = indices[s]; /**< actual subproblem index */

		/** cutoff = dual variable corresponding to the convex-combination constraint */
		double cutoff = dualsol_[sind];
		DSPdebugMessage("pricing out: %e < %e ? (colobj %e, status %d)\n", objs[s], cutoff, cxs[s], statuses[s]);

		/** retrieve subproblem solution */
		const CoinPackedVector* x = sols[s];

		/** create a column objective */
		double newcoef = cxs[s];

		/** take A x^k */
		mat_orig_->times(*x, Ax);

		/** clear a column vector */
		colvec.clear();

		/** convex combination constraints */
		if (statuses[s] != DSP_STAT_DUAL_INFEASIBLE)
			colvec.insert(sind, 1.0);

		/** original constraints */
		for (int i = 0; i < nrows_orig_; ++i)
			if (fabs(Ax[i]) > 1.0e-10)
				colvec.insert(nrows_conv_+i, Ax[i]);

		/** branching constraints */
#ifdef USE_ROW_TO_COL
		for (int i = 0; i < nrows_branch_; ++i) {
			int j = branch_row_to_col_[nrows_core_ + i];
			int sparse_index = x->findIndex(j);
			if (sparse_index == -1) continue;
			double val = x->getElements()[sparse_index];
			if (fabs(val) > 1.0e-10)
				colvec.insert(nrows_core_ + i, val);
		}
#else
		if (nrows_branch_ > 0) {
			double* dx = x->denseVector(ncols_orig_);
			for (int i = 0; i < nrows_branch_; ++i) {
				double val = branch_row_to_vec_[nrows_core_+i].dotProduct(dx);
				if (fabs(val) > 1.0e-10)
					colvec.insert(nrows_core_ + i, val);
			}
			FREE_ARRAY_PTR(dx);
		}
#endif

		if (statuses[s] == DSP_STAT_DUAL_INFEASIBLE || objs[s] < cutoff - 1.0e-4) {
			/** add the column vector */
			if (phase_ == 1)
				getSiPtr()->addCol(colvec.getNumElements(), colvec.getIndices(), colvec.getElements(), 
					0.0, COIN_DBL_MAX, 0.0);
			else if (phase_ == 2)
				getSiPtr()->addCol(colvec.getNumElements(), colvec.getIndices(), colvec.getElements(), 
					0.0, COIN_DBL_MAX, newcoef);

			/** store columns */
			cols_generated_.push_back(new DwCol(sind, getSiPtr()->getNumCols() - 1, *x, colvec, newcoef, 0.0, COIN_DBL_MAX));
			ngenerated_++;
		} else {
			/** store columns */
			cols_generated_.push_back(new DwCol(sind, -1, *x, colvec, newcoef, 0.0, COIN_DBL_MAX, false));
			DSPdebugMessage("added inactive column\n");
		}
	}
	DSPdebugMessage("Number of columns in the pool: %lu\n", cols_generated_.size());

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

DSP_RTN_CODE DwMaster::getLagrangianBound(
		std::vector<double>& objs /**< [in] subproblem objective values */) {
	BGN_TRY_CATCH

	/** calculate lower bound */
	dualobj_ = 0.0;
	const double* rlbd = getSiPtr()->getRowLower();
	const double* rubd = getSiPtr()->getRowUpper();
	for (int j = nrows_conv_; j < nrows_; ++j) {
		if (rlbd[j] > -1.0e+20)
			dualobj_ += dualsol_[j] * rlbd[j];
		else if (rubd[j] < 1.0e+20)
			dualobj_ += dualsol_[j] * rubd[j];
	}
	for (unsigned int s = 0; s < objs.size(); ++s) {
		dualobj_ += objs[s];
		DSPdebugMessage("subobj[%d] %e\n", s, objs[s]);
	}
	DSPdebugMessage("lb %e\n", dualobj_);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

bool DwMaster::terminationTest() {

	if (ngenerated_ == 0 || (!useBarrier_ && getSiPtr()->getIterationCount() == 0))
		return true;

	bool term = false;

	int status;
	if (phase_ == 1) {
		status = osi_->status();
		if (status == DSP_STAT_OPTIMAL && primobj_ < feastol_) {
			status_ = DSP_STAT_FEASIBLE;
			DSPdebugMessage("Phase 1 found a feasible solution! %e < %e\n", osi_->getPrimObjValue(), feastol_);
			term = true;
		}
	}
	if (phase_ == 2) {
		if (iterlim_ <= itercnt_) {
			status_ = DSP_STAT_LIM_ITERorTIME;
			term = true;
		} else if (bestdualobj_ >= bestprimobj_) {
			status_ = DSP_STAT_FEASIBLE;
			term = true;
		} else if (getRelApproxGap() < par_->getDblParam("DW/GAPTOL")) {
			status_ = DSP_STAT_OPTIMAL;
			term = true;
		}
	}

	return term;
}

DSP_RTN_CODE DwMaster::updateModel() {
	BGN_TRY_CATCH

	/** update the best dual objective */
	if (bestdualobj_ < dualobj_) bestdualobj_ = dualobj_;

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

DSP_RTN_CODE DwMaster::calculatePiA(
		std::vector<double>& piA /**< [out] pi^T A */) {
	BGN_TRY_CATCH

	piA.resize(ncols_orig_);
	mat_orig_->transposeTimes(&dualsol_[nrows_conv_], &piA[0]);
#ifdef USE_ROW_TO_COL
	for (int i = nrows_core_; i < nrows_; ++i)
		piA[branch_row_to_col_[i]] += dualsol_[i];
#else
	for (int i = 0; i < nrows_branch_; ++i) {
		const CoinPackedVector* vec = &branch_row_to_vec_[nrows_core_+i];
		int numels = vec->getNumElements();
		for (int j = 0; j < numels; ++j)
			piA[vec->getIndices()[j]] += vec->getElements()[j] * dualsol_[nrows_core_ + i];
	}
#endif

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

void DwMaster::setBranchingObjects(const DspBranchObj* branchobj) {
	/** shouldn't be null */
	if (branchobj == NULL)
		return;

	BGN_TRY_CATCH

	/** remove all columns */
	removeAllCols();

	/** remove all the branching rows */
	removeBranchingRows();

	/** add branching rows and columns */
	addBranchingRowsCols(branchobj);

	/** restore column bounds */
	clbd_node_ = clbd_orig_;
	cubd_node_ = cubd_orig_;

	/** update column bounds at the current node */
	int numColBranches = 0;
	for (unsigned j = 0; j < branchobj->getNumObjs(); ++j) {
		if (branchobj->getVector(j)->getNumElements() == 1) {
			//printf("Updated column bound [%d]: [%e, %e]\n", branchobj->getIndex(j), branchobj->getLb(j), branchobj->getUb(j));
			clbd_node_[branchobj->getIndex(j)] = branchobj->getLb(j);
			cubd_node_[branchobj->getIndex(j)] = branchobj->getUb(j);
			numColBranches++;
		}
	}

#ifdef DSP_DEBUG
	for (int j = 0; j < ncols_orig_; ++j) {
		if (fabs(clbd_node_[j] - clbd_orig_[j]) > 1.0e-8 || fabs(cubd_node_[j] - cubd_orig_[j]) > 1.0e-8)
			printf("Changed column bound [%d]: [%e,%e]\n", j, clbd_node_[j], cubd_node_[j]);
	}
#endif

	/** apply column bounds */
	if (numColBranches > 0) {
		std::vector<int> ncols_inds(ncols_orig_);
		CoinIotaN(&ncols_inds[0], ncols_orig_, 0);
		worker_->setColBounds(ncols_orig_, &ncols_inds[0], &clbd_node_[0], &cubd_node_[0]);
	}

	/** TODO: This is very specific to a certain SMIP form 
		and thus needs to be generalized later. */
	if (par_->getIntParam("DW/BRANCH") == BRANCH_DISJUNCTION_TEST) {
		if (model_->isStochastic()) {
			TssModel* tss = dynamic_cast<TssModel*>(model_);
			bool isFirstStage = true;
			for (unsigned j = 0; j < branchobj->getNumObjs(); ++j) {
				if (branchobj->getVector(j)->getMaxIndex() >= tss->getNumScenarios() * tss->getNumCols(0)) {
					isFirstStage = false;
					break;
				}
			}
			if (isFirstStage) {
				//printf("Adding general branching disjunctions (num: %d) to subproblem.\n", branchobj->getNumObjs());

				/** remove all rows added in the previous iteration */
				worker_->removeAddedRows();

				/** add general branching disjunctions */
				/** FIXME: add one by one is not efficient particularly in parallel mode. */
				CoinPackedVector* vec_copy = new CoinPackedVector;
				for (unsigned j = 0; j < branchobj->getNumObjs(); ++j) {
					if (branchobj->getVector(j)->getNumElements() > 1) {
						vec_copy->reserve(branchobj->getVector(j)->getNumElements());
						//printf("Reserved vec_copy %d\n", branchobj->getVector(j)->getNumElements());
						/** adjust column indices */
						for (int i = 0; i < branchobj->getVector(j)->getNumElements(); ++i)
							vec_copy->insert(
								branchobj->getVector(j)->getIndices()[i] % tss->getNumCols(0), 
								branchobj->getVector(j)->getElements()[i]);

						//printf("Branching object %d:\n", j);
						//DspMessage::printArray(vec_copy);
						//printf("lb %e ub %e\n", branchobj->getLb(j), branchobj->getUb(j));
						/** add branching disjunction */
						worker_->addRow(vec_copy, branchobj->getLb(j), branchobj->getUb(j));
						vec_copy->clear();
					}
				}
				FREE_PTR(vec_copy);
			} else {
				message_->print(0, "General disjunction branching is allowed on first-stage variables only.\n");
			}
		} else {
			message_->print(0, "General disjunction branching is allowed on SMIP only.\n");
		}
	}

	/** set known best bound */
	bestdualobj_ = COIN_DBL_MAX;

	/** set braching object pointer */
	branchObj_ = branchobj;

	END_TRY_CATCH(;)
}

void DwMaster::removeAllCols() {
	int ndelcols = getSiPtr()->getNumCols();
	std::vector<int> delcols;
	delcols.reserve(ndelcols);
	for (int j = 0; j < ndelcols; ++j)
		delcols.push_back(j);
	getSiPtr()->deleteCols(ndelcols, &delcols[0]);
	DSPdebugMessage("Deleted %d columns in the master.\n", ndelcols);
}

void DwMaster::removeBranchingRows() {
	if (nrows_branch_ > 0) {
		std::vector<int> delrows;
		delrows.reserve(nrows_branch_);
		for (int i = nrows_core_; i < nrows_; ++i)
			delrows.push_back(i);
		getSiPtr()->deleteRows(nrows_branch_, &delrows[0]);
		DSPdebugMessage("Deleted %d rows in the master.\n", nrows_branch_);
		nrows_branch_ = 0;
	}
}

void DwMaster::addBranchingRowsCols(const DspBranchObj* branchobj) {

	/** count nrows_branch_ */
	if (par_->getBoolParam("DW/MASTER/BRANCH_ROWS")) {
		for (int j = 0; j < branchobj->getNumObjs(); ++j) {
#ifdef USE_ROW_TO_COL
			if (branchobj->getLb(j) > clbd_orig_[branchobj->getIndex(j)]) {
				branch_row_to_col_[nrows_core_ + nrows_branch_] = branchobj->getIndex(j);
				//getSiPtr()->addRow(0, NULL, NULL, branchobj->getLb(j), COIN_DBL_MAX);
				addBranchingRow(branchobj->getLb(j), COIN_DBL_MAX);
				nrows_branch_++;
			}
			if (branchobj->getUb(j) < cubd_orig_[branchobj->getIndex(j)]) {
				branch_row_to_col_[nrows_core_ + nrows_branch_] = branchobj->getIndex(j);
				//getSiPtr()->addRow(0, NULL, NULL, -COIN_DBL_MAX, branchobj->getUb(j));
				addBranchingRow(-COIN_DBL_MAX, branchobj->getUb(j));
				nrows_branch_++;
			}
#else
			if (branchobj->getLb(j) > clbd_orig_[branchobj->getIndex(j)]) {
				branch_row_to_vec_[nrows_core_ + nrows_branch_] = *(branchobj->getVector(j));
				addBranchingRow(branchobj->getLb(j), COIN_DBL_MAX);
				nrows_branch_++;
			}
			if (branchobj->getUb(j) < cubd_orig_[branchobj->getIndex(j)]) {
				branch_row_to_vec_[nrows_core_ + nrows_branch_] = *(branchobj->getVector(j));
				addBranchingRow(-COIN_DBL_MAX, branchobj->getUb(j));
				nrows_branch_++;
			}
#endif
		}
	}

	/** update number of rows */
	nrows_ = nrows_core_ + nrows_branch_;
	DSPdebugMessage("nrows_ %d nrows_core_ %d nrows_branch_ %d\n", nrows_, nrows_core_, nrows_branch_);

	if (par_->getBoolParam("DW/MASTER/REUSE_COLS")) {
		/** adding columns */
		std::vector<int> col_inds;
		std::vector<double> col_elems;

		/** add branching rows */
		for (auto it = cols_generated_.begin(); it != cols_generated_.end(); it++) {
			(*it)->active_ = true;
			(*it)->age_ = 0;
			double* dx = (*it)->x_.denseVector(ncols_orig_);
			for (unsigned j = 0; j < branchobj->getNumObjs(); ++j) {
#ifdef USE_ROW_TO_COL
				int sparse_index = (*it)->x_.findIndex(branchobj->getIndex(j));
				double val = 0.0;
				if (sparse_index <= -1) continue;
				val = (*it)->x_.getElements()[sparse_index];
#else
				double val = 0.0;
				val = branch_row_to_vec_[nrows_core_+j].dotProduct(dx);
#endif
				if (val < branchobj->getLb(j) || val > branchobj->getUb(j)) {
					(*it)->active_ = false;
					(*it)->age_ = COIN_INT_MAX;
					break;
				}
			}
			FREE_ARRAY_PTR(dx);

			if ((*it)->active_) {
				if (par_->getBoolParam("DW/MASTER/BRANCH_ROWS")) {
					DSP_RTN_CHECK_THROW(updateCol(*it));
				} else {
					/** create a column for core rows */
					col_inds.clear();
					col_elems.clear();
					col_inds.reserve((*it)->col_.getNumElements());
					col_elems.reserve((*it)->col_.getNumElements());
					for (int i = 0; i < (*it)->col_.getNumElements(); ++i) {
						if ((*it)->col_.getIndices()[i] < nrows_core_) {
							col_inds.push_back((*it)->col_.getIndices()[i]);
							col_elems.push_back((*it)->col_.getElements()[i]);
						}
					}

					/** assign the core-row column */
					(*it)->col_.setVector(col_inds.size(), &col_inds[0], &col_elems[0]);
				}
				/** set master index */
				(*it)->master_index_ = getSiPtr()->getNumCols();

				/** add column */
				//getSiPtr()->addCol((*it)->col_, 0.0, COIN_DBL_MAX, (*it)->obj_);
				addBranchingCol((*it)->col_, (*it)->obj_);
			}
		}
		message_->print(2, "Appended dynamic columns in the master (%d / %u cols).\n", getSiPtr()->getNumRows(), cols_generated_.size());
	} else {
		for (unsigned i = 0; i < cols_generated_.size(); ++i)
			FREE_PTR(cols_generated_[i]);
		cols_generated_.clear();
	}
}

void DwMaster::addBranchingRow(double lb, double ub) {
	getSiPtr()->addRow(0, NULL, NULL, lb, ub);
}

void DwMaster::addBranchingCol(const CoinPackedVector& col, double obj) {
	getSiPtr()->addCol(col.getNumElements(), col.getIndices(), col.getElements(), 
		0.0, COIN_DBL_MAX, obj);
}

void DwMaster::printIterInfo() {
	message_->print(2, "[Phase %d] Iteration %3d: Master objective %e, ", phase_, itercnt_, primobj_);
	if (phase_ == 2) {
		if (bestdualobj_ > -1.0e+50)
			message_->print(2, "Lb %e (gap %.2f %%), ", bestdualobj_, getRelApproxGap()*100);
		else
			message_->print(2, "Lb -Inf, ");
	}
	message_->print(2, "nrows %d, ncols %d, ", getSiPtr()->getNumRows(), getSiPtr()->getNumCols());
	if (!useBarrier_)
		message_->print(2, "itercnt %d, ", getSiPtr()->getIterationCount());
	message_->print(2, "timing (total %.2f, master %.2f, gencols %.2f), statue %d\n",
			t_total_, t_master_, t_colgen_, status_);
}

DSP_RTN_CODE DwMaster::switchToPhase2() {
	BGN_TRY_CATCH
	if (phase_ == 1) {
		/** delete auxiliary columns */
		getSiPtr()->deleteCols(auxcolindices_.size(), &auxcolindices_[0]);
		auxcolindices_.clear();
		DSPdebugMessage("Phase 2 has %d rows and %d columns.\n", getSiPtr()->getNumRows(), getSiPtr()->getNumCols());

		/** set objective function coefficients */
		for (unsigned k = 0, j = 0; k < cols_generated_.size(); ++k)
			if (cols_generated_[k]->active_) {
				if (j >= getSiPtr()->getNumCols()) {
					message_->print(0, "Trying to access invalid column index %d (ncols %d)\n", j, getSiPtr()->getNumCols());
					return DSP_RTN_ERR;
				}
				getSiPtr()->setObjCoeff(j, cols_generated_[k]->obj_);
				j++;
			}
		phase_ = 2;
	}
	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

DspOsi * DwMaster::createDspOsi() {
	DspOsi* osi = NULL;

	BGN_TRY_CATCH

	switch(par_->getIntParam("DW/MASTER/SOLVER")) {
		case OsiCpx:
#ifdef DSP_HAS_CPX
			osi = new DspOsiCpx();
#else
			throw CoinError("Cplex is not available.", "createDspOsi", "DwMaster");
#endif
			break;
		case OsiGrb:
#ifdef DSP_HAS_GRB
			osi = new DspOsiGrb();
#else
			throw CoinError("Gurobi is not available.", "createDspOsi", "DwMaster");
#endif
			break;
		case OsiClp:
			osi = new DspOsiClp();
			break;
		default:
			throw CoinError("Invalid parameter value", "initDualSolver", "DwMaster");
			break;
	}

	END_TRY_CATCH_RTN(;,osi)
	
	return osi;
}
