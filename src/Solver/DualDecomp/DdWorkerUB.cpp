/*
 * DdWorkerUB.cpp
 *
 *  Created on: Mar 28, 2016
 *      Author: kibaekkim
 */

//#define DSP_DEBUG
#include "Model/TssModel.h"
#include <Solver/DualDecomp/DdWorkerUB.h>
#include "SolverInterface/SolverInterfaceClp.h"

#ifndef NO_CPX
	#include "SolverInterface/SolverInterfaceCpx.h"
#endif

#ifndef NO_SCIP
	#include "SolverInterface/SolverInterfaceScip.h"
	#include "Solver/DualDecomp/SCIPconshdlrBendersDd.h"
	#include "SolverInterface/SCIPbranchruleLB.h"
#endif

DdWorkerUB::DdWorkerUB(
		DspParams * par,
		DecModel * model,
		DspMessage * message):
DdWorker(par, model, message),
mat_mp_(NULL),
rlbd_org_(NULL),
rubd_org_(NULL),
si_(NULL),
objvals_(NULL),
statuses_(NULL),
ub_(0.0) {}

DdWorkerUB::~DdWorkerUB() {
	FREE_2D_PTR(par_->getIntPtrParamSize("ARR_PROC_IDX"), mat_mp_);
	FREE_2D_PTR(par_->getIntPtrParamSize("ARR_PROC_IDX"), rlbd_org_);
	FREE_2D_PTR(par_->getIntPtrParamSize("ARR_PROC_IDX"), rubd_org_);
	FREE_2D_PTR(par_->getIntPtrParamSize("ARR_PROC_IDX"), si_);
	FREE_ARRAY_PTR(objvals_);
	FREE_ARRAY_PTR(statuses_);
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
	FREE_ARRAY_PTR(obj_reco)

	/** recourse problem data */
	CoinPackedMatrix * mat_reco = NULL;
	double * clbd_reco   = NULL;
	double * cubd_reco   = NULL;
	double * obj_reco    = NULL;
	char *   ctype_reco  = NULL;

	TssModel* tss = NULL;

	BGN_TRY_CATCH

	int nsubprobs = par_->getIntPtrParamSize("ARR_PROC_IDX");
	tss = dynamic_cast<TssModel*>(model_);
	if (tss == NULL)
		throw "This is not a stochastic programming problem.";

	/** allocate memory */
	mat_mp_    = new CoinPackedMatrix * [nsubprobs];
	rlbd_org_  = new double * [nsubprobs];
	rubd_org_  = new double * [nsubprobs];
	si_        = new SolverInterface * [nsubprobs];
	objvals_   = new double [nsubprobs];
	statuses_  = new int [nsubprobs];

	for (int s = 0; s < nsubprobs; ++s) {

		/** copy recourse problem */
		DSP_RTN_CHECK_THROW(tss->copyRecoProb(par_->getIntPtrParam("ARR_PROC_IDX")[s],
				mat_mp_[s], mat_reco, clbd_reco, cubd_reco, ctype_reco,
				obj_reco, rlbd_org_[s], rubd_org_[s]));

		/** creating solver interface */
    	switch (par_->getIntParam("SOLVER/MIP")) {
    	case CPLEX:
#ifndef NO_CPX
    		si_[s] = new SolverInterfaceCpx(par_);
    		break;
#endif
    	case SCIP:
#ifndef NO_SCIP
            si_[s] = new SolverInterfaceScip(par_);
            break;
#endif
    	default:
    		break;
    	}

	    /** no display */
	    si_[s]->setPrintLevel(0);

	    /** load problem */
	    si_[s]->loadProblem(mat_reco, clbd_reco, cubd_reco, obj_reco, ctype_reco, rlbd_org_[s], rubd_org_[s]);
	    DSPdebug(mat_reco->verifyMtx(4));
    }
	END_TRY_CATCH_RTN(FREE_MEMORY, DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
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

	int nrows = mat_mp_[0]->getNumRows(); /** retrieve the number of rows in subproblem */
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
		const double* rlbd = si_[s]->getRowLower();
		const double* rubd = si_[s]->getRowUpper();
		for (int i = si_[s]->getNumRows() - 1; i >= 0; --i) {
			if (rlbd_org_[s][i] > -COIN_DBL_MAX)
				si_[s]->setRowLower(i, rlbd[i] - Tx[i]);
			if (rubd_org_[s][i] < COIN_DBL_MAX)
				si_[s]->setRowUpper(i, rubd[i] - Tx[i]);
		}

		/** first-stage objective value */
		cx_weighted += cx * tss->getProbability()[par_->getIntPtrParam("ARR_PROC_IDX")[s]];

		FREE_ARRAY_PTR(Tx)
	}

	DSP_RTN_CHECK_RTN_CODE(solve());
	ub_ += cx_weighted;

	/** restore row bounds */
	for (int s = nsubprobs - 1; s >= 0; --s) {
		for (int i = si_[s]->getNumRows() - 1; i >= 0; --i) {
			si_[s]->setRowLower(i, rlbd_org_[s][i]);
			si_[s]->setRowUpper(i, rubd_org_[s][i]);
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
		si_[s]->setTimeLimit(
				CoinMin(CoinMax(0.01, time_remains_),
				par_->getDblParam("MIP/TIME_LIM")));

		/** solve */
		si_[s]->solve();

		/** check status. there might be unexpected results. */
		switch (si_[s]->getStatus()) {
		case DSP_STAT_OPTIMAL:
		case DSP_STAT_LIM_ITERorTIME:
		case DSP_STAT_STOPPED_GAP:
		case DSP_STAT_STOPPED_NODE:
		case DSP_STAT_STOPPED_TIME:
			break;
		default:
			status_ = DSP_STAT_MW_STOP;
			message_->print(0,
					"Warning: subproblem %d solution status is %d\n", s,
					si_[s]->getStatus());
			break;
		}
		if (status_ == DSP_STAT_MW_STOP)
			break;

		primobj += si_[s]->getPrimalBound();
		dualobj += si_[s]->getDualBound();
		total_cputime += CoinCpuTime() - cputime;
		total_walltime += CoinGetTimeOfDay() - walltime;

		/** consume time */
		time_remains_ -= CoinGetTimeOfDay() - walltime;
	}

	/** get primal objective */
	ub_ = primobj;

	/** update statistics */
	s_statuses_.push_back(status_);
	s_primobjs_.push_back(primobj);
	s_dualobjs_.push_back(dualobj);
	s_cputimes_.push_back(total_cputime);
	s_walltimes_.push_back(total_walltime);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}
