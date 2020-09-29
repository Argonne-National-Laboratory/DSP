/*
 * DdWorkerUB.cpp
 *
 *  Created on: Mar 28, 2016
 *      Author: kibaekkim
 */

// #define DSP_DEBUG
#include "Model/TssModel.h"
#include "Solver/DualDecomp/DdWorkerUB.h"
#include "SolverInterface/DspOsiCpx.h"
#include "SolverInterface/DspOsiGrb.h"
#include "SolverInterface/DspOsiScip.h"

#ifdef DSP_HAS_SCIP
#include "Solver/DualDecomp/SCIPconshdlrBendersDd.h"
#endif

DdWorkerUB::DdWorkerUB(
		DecModel *   model,  /**< model pointer */
		DspParams *  par,    /**< parameter pointer */
		DspMessage * message /**< message pointer */):
DdWorker(model, par, message),
bestub_(COIN_DBL_MAX),
mat_mp_(NULL),
rlbd_org_(NULL),
rubd_org_(NULL),
osi_(NULL),
osi_dro_(NULL),
ub_(0.0) {}

DdWorkerUB::DdWorkerUB(const DdWorkerUB& rhs) : 
DdWorker(rhs),
bestub_(rhs.bestub_),
primsols_(rhs.primsols_),
ub_(rhs.ub_) {
	// number of subproblems to take care of
	int nsubprobs = par_->getIntPtrParamSize("ARR_PROC_IDX");

	// allocate memory
	mat_mp_    = new CoinPackedMatrix * [nsubprobs];
	rlbd_org_  = new double * [nsubprobs];
	rubd_org_  = new double * [nsubprobs];
	osi_        = new DspOsi * [nsubprobs];
	for (int s = 0; s < nsubprobs; ++s) {
		mat_mp_[s] = new CoinPackedMatrix(*(rhs.mat_mp_[s]));
		rlbd_org_[s] = new double [osi_[s]->si_->getNumRows()];
		rubd_org_[s] = new double [osi_[s]->si_->getNumRows()];
		osi_[s] = rhs.osi_[s]->clone();
		CoinCopyN(rhs.rlbd_org_[s], osi_[s]->si_->getNumRows(), rlbd_org_[s]);
		CoinCopyN(rhs.rubd_org_[s], osi_[s]->si_->getNumRows(), rubd_org_[s]);
	}
	if (model_->isDro())
		osi_dro_ = rhs.osi_dro_->clone();
}

DdWorkerUB::~DdWorkerUB() {
	FREE_2D_PTR(par_->getIntPtrParamSize("ARR_PROC_IDX"), mat_mp_);
	FREE_2D_ARRAY_PTR(par_->getIntPtrParamSize("ARR_PROC_IDX"), rlbd_org_);
	FREE_2D_ARRAY_PTR(par_->getIntPtrParamSize("ARR_PROC_IDX"), rubd_org_);
	FREE_2D_PTR(par_->getIntPtrParamSize("ARR_PROC_IDX"), osi_);
	FREE_PTR(osi_dro_);
}

DSP_RTN_CODE DdWorkerUB::init() {
	BGN_TRY_CATCH
	/** status */
	status_ = DSP_STAT_MW_CONTINUE;
	/** create problem */
	DSP_RTN_CHECK_THROW(createProblem());
	END_TRY_CATCH_RTN(;, DSP_RTN_ERR)
	return DSP_RTN_OK;
}

DSP_RTN_CODE DdWorkerUB::createProblem() {
#define FREE_MEMORY            \
	FREE_PTR(mat_reco)         \
	FREE_ARRAY_PTR(clbd_reco)  \
	FREE_ARRAY_PTR(cubd_reco)  \
	FREE_ARRAY_PTR(ctype_reco) \
	FREE_ARRAY_PTR(obj_reco)   \
	FREE_ARRAY_PTR(clbd_dro)   \
	FREE_ARRAY_PTR(cubd_dro)   \
	FREE_ARRAY_PTR(rlbd_dro)   \
	FREE_ARRAY_PTR(rubd_dro)   \
	FREE_ARRAY_PTR(obj_dro)    \
	FREE_ARRAY_PTR(elem_dro)   \
	FREE_ARRAY_PTR(ind_dro)    \
	FREE_ARRAY_PTR(bgn_dro)    \
	FREE_ARRAY_PTR(len_dro)

	/** recourse problem data */
	CoinPackedMatrix * mat_reco = NULL;
	double * clbd_reco   = NULL;
	double * cubd_reco   = NULL;
	double * obj_reco    = NULL;
	char *   ctype_reco  = NULL;

	/** DRO UB problem */
	double * clbd_dro = NULL;
	double * cubd_dro = NULL;
	double * obj_dro = NULL;
	double * rlbd_dro = NULL;
	double * rubd_dro = NULL;
	double * elem_dro = NULL;
	int * ind_dro = NULL;
	int * bgn_dro = NULL;
	int * len_dro = NULL;

	TssModel* tss = NULL;

	bool has_integer = false;

	BGN_TRY_CATCH

	int nsubprobs = par_->getIntPtrParamSize("ARR_PROC_IDX");
	tss = dynamic_cast<TssModel*>(model_);
	if (tss == NULL)
		throw "This is not a stochastic programming problem.";

	/** allocate memory */
	mat_mp_    = new CoinPackedMatrix * [nsubprobs];
	rlbd_org_  = new double * [nsubprobs];
	rubd_org_  = new double * [nsubprobs];
	osi_       = new DspOsi * [nsubprobs];
	primsols_.resize(nsubprobs);

	for (int s = 0; s < nsubprobs; ++s) {

		/** copy recourse problem */
		DSP_RTN_CHECK_THROW(model_->copyRecoProb(par_->getIntPtrParam("ARR_PROC_IDX")[s],
				mat_mp_[s], mat_reco, clbd_reco, cubd_reco, ctype_reco,
				obj_reco, rlbd_org_[s], rubd_org_[s]));

		if (model_->isDro()) {
			for (int j = 0; j < tss->getNumCols(1); ++j) {
				obj_reco[j] /= tss->getProbability()[par_->getIntPtrParam("ARR_PROC_IDX")[s]];
			}
		}

		for (int j = 0; j < mat_reco->getNumCols(); ++j)
			if (ctype_reco[j] != 'C') {
				has_integer = true;
				break;
			}

		/** creating solver interface */
		osi_[s] = createDspOsi();
		if (!osi_[s]) throw CoinError("Failed to create DspOsi", "createProblem", "DdWorkerUB");

	    /** no display */
	    osi_[s]->setLogLevel(0);
		DSPdebug(osi_[s]->setLogLevel(5));

	    /** load problem */
	    osi_[s]->si_->loadProblem(*mat_reco, clbd_reco, cubd_reco, obj_reco, rlbd_org_[s], rubd_org_[s]);
		for (int j = 0; j < mat_reco->getNumCols(); ++j) {
			if (ctype_reco[j] != 'C')
				osi_[s]->si_->setInteger(j);
		}
	    DSPdebug(mat_reco->verifyMtx(4));
		DSPdebugMessage("number of integers: %d\n", osi_[s]->si_->getNumIntegers());

		/** allocate array size for each scenario primal solution */
		primsols_[s].resize(osi_[s]->si_->getNumCols());

		FREE_MEMORY
    }

	/** create DRO upper bounding problem */
	if (model_->isDro()) {
		int ncols_dro = 1 + model_->getNumReferences();
		int nrows_dro = model_->getNumReferences() * tss->getNumScenarios();
		int nzcnt_dro = nrows_dro * 2;

		/** allocate memory */
		clbd_dro = new double [ncols_dro];
		cubd_dro = new double [ncols_dro];
		obj_dro = new double [ncols_dro];
		rlbd_dro = new double [nrows_dro];
		rubd_dro = new double [nrows_dro];
		bgn_dro = new int [nrows_dro+1];
		len_dro = new int [nrows_dro];
		ind_dro = new int [nzcnt_dro];
		elem_dro = new double [nzcnt_dro];

		/** column bounds */
		clbd_dro[0] = 0.0;
		CoinFillN(clbd_dro + 1, ncols_dro - 1, -COIN_DBL_MAX);
		CoinFillN(cubd_dro, ncols_dro, COIN_DBL_MAX);

		/** row bounds */
		CoinFillN(rlbd_dro, nrows_dro, -COIN_DBL_MAX);
		CoinFillN(rubd_dro, nrows_dro, COIN_DBL_MAX);

		/** objective coefficients */
		obj_dro[0] = model_->getWassersteinSize();
		for (int j = 1; j < ncols_dro; ++j)
			obj_dro[j] = model_->getReferenceProbability(j-1);

		int pos_dro = 0;
		int rnum = 0;
		for (int k = 0; k < tss->getNumScenarios(); ++k) {
			for (int s = 0; s < model_->getNumReferences(); ++s) {
				rnum = k * model_->getNumReferences() + s;

				bgn_dro[rnum] = pos_dro;

				// alpha
				ind_dro[pos_dro] = 0; 
				elem_dro[pos_dro] = model_->getWassersteinDist(s,k);
				pos_dro++;

				// beta_s
				ind_dro[pos_dro] = 1 + s; 
				elem_dro[pos_dro] = 1.0;
				pos_dro++;

				len_dro[rnum] = pos_dro - bgn_dro[rnum];
			}
		}
		bgn_dro[nrows_dro] = pos_dro;
		assert(pos_dro == nzcnt_dro);

		/** constraint matrix */
		CoinPackedMatrix *mat_dro = new CoinPackedMatrix(false, ncols_dro, nrows_dro, nzcnt_dro, elem_dro, ind_dro, bgn_dro, len_dro);
		DSPdebug(mat_dro->verifyMtx(4));
		
		/** creating solver interface */
		osi_dro_ = createDspOsi();
		if (!osi_dro_) throw CoinError("Failed to create DspOsi", "createProblem", "DdWorkerUB");

	    /** no display */
	    osi_dro_->setLogLevel(0);

	    /** load problem */
	    osi_dro_->si_->loadProblem(*mat_dro, clbd_dro, cubd_dro, obj_dro, rlbd_dro, rubd_dro);

		/** free matrix */
		FREE_PTR(mat_dro);
	}

	END_TRY_CATCH_RTN(FREE_MEMORY, DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

double DdWorkerUB::evaluate(int n, double* solution) {
	std::vector<int> indices;
	std::vector<double> elements;
	for (int i = 0; i < n; ++i)
		if (fabs(solution[i]) > 1.0e-10) {
			indices.push_back(i);
			elements.push_back(solution[i]);
		}

	CoinPackedVector *s = new CoinPackedVector(indices.size(), &indices[0], &elements[0]);
	double ub = evaluate(s);
	
	delete s;

	return ub;
}

double DdWorkerUB::evaluate(CoinPackedVector* solution) {
#define FREE_MEMORY \
	FREE_ARRAY_PTR(Tx) \
	FREE_ARRAY_PTR(x)

	double * Tx = NULL;
	double* x = NULL;
	TssModel* tss = NULL;

	BGN_TRY_CATCH

	tss = dynamic_cast<TssModel*>(model_);
	if (tss == NULL)
		throw "This is not a stochastic programming problem.";

	int nsubprobs = par_->getIntPtrParamSize("ARR_PROC_IDX");

	/** first-stage objective value */
	double cx = solution->dotProduct(tss->getObjCore(0));
	double cx_weighted = 0.0;

	/** allocate memory */
	x = solution->denseVector(tss->getNumCols(0));
	for (int s = nsubprobs - 1; s >= 0; --s) {
		/** calculate Tx */
		Tx = new double [mat_mp_[s]->getNumRows()];
		mat_mp_[s]->times(x, Tx);

		/** adjust row bounds */
		const double* rlbd = osi_[s]->si_->getRowLower();
		const double* rubd = osi_[s]->si_->getRowUpper();
		for (int i = osi_[s]->si_->getNumRows() - 1; i >= 0; --i) {
			if (rlbd_org_[s][i] > -COIN_DBL_MAX)
				osi_[s]->si_->setRowLower(i, rlbd[i] - Tx[i]);
			if (rubd_org_[s][i] < COIN_DBL_MAX)
				osi_[s]->si_->setRowUpper(i, rubd[i] - Tx[i]);
		}

		/** first-stage objective value */
		cx_weighted += cx * tss->getProbability()[par_->getIntPtrParam("ARR_PROC_IDX")[s]];

		FREE_ARRAY_PTR(Tx)
	}
	DSPdebugMessage("cx_weighted %e\n", cx_weighted);

	/** set initial status */
	status_ = DSP_STAT_MW_CONTINUE;

	/** solve upper bounding problem */
	DSP_RTN_CHECK_RTN_CODE(solve());

	/** update upper bound */
	if (status_ == DSP_STAT_MW_CONTINUE)
		ub_ += cx_weighted;
	else
		ub_ = COIN_DBL_MAX;

	/** restore row bounds */
	for (int s = nsubprobs - 1; s >= 0; --s) {
		for (int i = osi_[s]->si_->getNumRows() - 1; i >= 0; --i) {
			osi_[s]->si_->setRowLower(i, rlbd_org_[s][i]);
			osi_[s]->si_->setRowUpper(i, rubd_org_[s][i]);
		}
	}

	END_TRY_CATCH_RTN(FREE_MEMORY,COIN_DBL_MAX)
	FREE_MEMORY

	return ub_;
#undef FREE_MEMORY
}

DSP_RTN_CODE DdWorkerUB::solve() {
	double cputime;
	double walltime;

	BGN_TRY_CATCH

	double primobj = 0.0;
	double dualobj = 0.0;
	double total_cputime = 0.0;
	double total_walltime = 0.0;
	int nsubprobs = par_->getIntPtrParamSize("ARR_PROC_IDX");

	for (unsigned s = 0; s < nsubprobs; ++s)
	{
		cputime = CoinCpuTime();
		walltime = CoinGetTimeOfDay();

		/** set time limit */
		osi_[s]->setTimeLimit(
				CoinMin(CoinMax(0.01, time_remains_),
				par_->getDblParam("DD/SUB/TIME_LIM")));

		/** solve */
		osi_[s]->solve();

		/** check status. there might be unexpected results. */
		int status = osi_[s]->status();
		DSPdebugMessage("status = %d\n", status);
		switch (status) {
		case DSP_STAT_OPTIMAL:
		case DSP_STAT_LIM_ITERorTIME:
		case DSP_STAT_STOPPED_GAP:
		case DSP_STAT_STOPPED_NODE:
		case DSP_STAT_STOPPED_TIME:
			break;
		default:
			status_ = DSP_STAT_MW_STOP;
			message_->print(10,
					"Warning: subproblem %d solution status is %d\n", s,
					status);
			break;
		}
		if (status_ == DSP_STAT_MW_STOP) {
			primobj = COIN_DBL_MAX;
			dualobj = -COIN_DBL_MAX;
			break;
		}

		primobj += osi_[s]->getPrimObjValue();
		dualobj += osi_[s]->getDualObjValue();
		CoinCopyN(osi_[s]->si_->getColSolution(), osi_[s]->si_->getNumCols(), &primsols_[s][0]);
		total_cputime += CoinCpuTime() - cputime;
		total_walltime += CoinGetTimeOfDay() - walltime;

		/** consume time */
		time_remains_ -= CoinGetTimeOfDay() - walltime;
	}

	/** solve the DRO UB problem */
	if (model_->isDro() && status_ != DSP_STAT_MW_STOP) {
		assert(nsubprobs == model_->getNumSubproblems());
		for (int k = 0; k < nsubprobs; ++k) {
			for (int s = 0; s < model_->getNumReferences(); ++s) {
				osi_dro_->si_->setRowLower(k * model_->getNumReferences() + s, osi_[k]->si_->getObjValue());
			}
		}
		osi_dro_->solve();

		if (osi_dro_->si_->isProvenOptimal()) {
			primobj = osi_dro_->getPrimObjValue();
		} else {
			primobj = COIN_DBL_MAX;
		}
	}

	/** get primal objective */
	ub_ = primobj;
	DSPdebugMessage("ub_ = %e\n", ub_);
	DSPdebugMessage("status_ %d\n", status_);

	/** update statistics */
	s_statuses_.push_back(status_);
	s_primobjs_.push_back(primobj);
	s_dualobjs_.push_back(dualobj);
	s_cputimes_.push_back(total_cputime);
	s_walltimes_.push_back(total_walltime);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DspOsi * DdWorkerUB::createDspOsi() {
	DspOsi * osi = NULL;
	BGN_TRY_CATCH

	switch (par_->getIntParam("DW/SUB/SOLVER")) {
	case OsiCpx:
#ifdef DSP_HAS_CPX
		osi = new DspOsiCpx();
		CPXsetintparam(dynamic_cast<DspOsiCpx*>(osi)->cpx_->getEnvironmentPtr(), CPX_PARAM_SCRIND, CPX_OFF);
#else
		throw CoinError("Cplex is not available.", "createDspOsi", "DdWorkerUB");
#endif
		break;
	case OsiGrb:
#ifdef DSP_HAS_GRB
		osi = new DspOsiGrb();
#else
		throw CoinError("Gurobi is not available.", "createDspOsi", "DdWorkerUB");
#endif
		break;
	case OsiScip:
#ifdef DSP_HAS_SCIP
		osi = new DspOsiScip();
#else
		throw CoinError("Scip is not available.", "createDspOsi", "DdWorkerUB");
#endif
		break;
	default:
		throw CoinError("Invalid parameter value", "createDspOsi", "DdWorkerUB");
		break;
	}

	END_TRY_CATCH(;)
	return osi;
}
