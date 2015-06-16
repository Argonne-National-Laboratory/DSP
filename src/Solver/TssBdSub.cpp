/*
 * TssBdSub.cpp
 *
 *  Created on: Dec 5, 2014
 *      Author: kibaekkim
 */

//#define DSP_DEBUG

#include "Utility/StoMessage.h"
#include "Solver/TssBdSub.h"
#include "Solver/SolverInterfaceSpx.h"
#include "Solver/SolverInterfaceClp.h"
#include "Solver/SolverInterfaceScip.h"

/** Coin */
#include "OsiSpxSolverInterface.hpp"
#include "OsiClpSolverInterface.hpp"

#define COIN_OSI OsiSpxSolverInterface
#define DSP_SI   SolverInterfaceSpx
#define SOLVER_FEASIBILITY_FIRST

//#define TSS_BD_SUB_DEBUG
//#define TSS_BD_SUB_DEBUG_MORE

TssBdSub::TssBdSub(StoParam * par):
	par_(par),
	nSubs_(0),
	nAuxvars_(1),
	excludedScenarios_(NULL),
	weight_(NULL),
	mat_mp_(NULL),
	recourse_(NULL),
	cglp_(NULL),
	warm_start_(NULL),
	cut_generation_wall_time_(0.0),
	nFeasSolutions_(0),
	status_(NULL)
{
	/** nothing to do */
}

TssBdSub::~TssBdSub()
{
	FREE_ARRAY_PTR(status_);
	FREE_ARRAY_PTR(excludedScenarios_);
	FREE_2D_PTR(nSubs_,mat_mp_);
	FREE_2D_PTR(nSubs_,recourse_);
	FREE_2D_PTR(nSubs_,warm_start_);
	FREE_2D_PTR(nSubs_,cglp_);
	nSubs_    = 0;
	nAuxvars_ = 1;
}

/** load problem */
STO_RTN_CODE TssBdSub::loadProblem(
		TssModel * model,
		int naugs,
		int * augs,
		int nAuxvars)
{
#define FREE_MEMORY            \
	FREE_PTR(mat_reco)         \
	FREE_ARRAY_PTR(clbd_reco)  \
	FREE_ARRAY_PTR(cubd_reco)  \
	FREE_ARRAY_PTR(ctype_reco) \
	FREE_ARRAY_PTR(obj_reco)   \
	FREE_ARRAY_PTR(rlbd_reco)  \
	FREE_ARRAY_PTR(rubd_reco)

	/** recourse problem data */
	CoinPackedMatrix * mat_reco = NULL;
	double * clbd_reco   = NULL;
	double * cubd_reco   = NULL;
	double * obj_reco    = NULL;
	char *   ctype_reco  = NULL;
	double * rlbd_reco   = NULL;
	double * rubd_reco   = NULL;
	assert(model);

	BGN_TRY_CATCH

	nSubs_    = model->getNumScenarios();
	weight_   = model->getProbability();
	mat_mp_   = new CoinPackedMatrix * [nSubs_];
	if (par_->TssDdCacheRecourse_)
		recourse_ = new SolverInterface * [nSubs_];
	cglp_     = new OsiSolverInterface * [nSubs_];
	status_   = new int [nSubs_];
	nAuxvars_ = nAuxvars;
	for (int s = nSubs_ - 1; s >= 0; --s)
	{
		mat_mp_[s] = NULL;
		if (par_->TssDdCacheRecourse_) recourse_[s] = NULL;
		cglp_[s] = NULL;
	}

#ifndef NO_WARMSTART_TEST
	/** initialize warm start information */
	warm_start_ = new CoinWarmStart * [nSubs_];
	for (int s = nSubs_ - 1; s >= 0; --s)
		warm_start_[s] = NULL;
#endif

	/** mark excluded scenarios */
	excludedScenarios_ = new bool [nSubs_];
	CoinFillN(excludedScenarios_, nSubs_, false);
	for (int s = 0; s < naugs; ++s)
		excludedScenarios_[augs[s]] = true;

	/** for included scenarios */
#if 0
	for (int s = 0; s < nSubs_; ++s)
	{
#else
	for (unsigned int i = 0; i < scenarios_.size(); ++i)
	{
		int s = scenarios_[i];
#endif
		if (excludedScenarios_[s]) continue;

		/** copy recourse problem */
		STO_RTN_CHECK_THROW(model->copyRecoProb(s, mat_mp_[s],
				mat_reco, clbd_reco, cubd_reco, ctype_reco,
				obj_reco, rlbd_reco, rubd_reco), "copyRecoProb", "TssModel");

		/** create cache for getting recourse solution */
		if (par_->TssDdCacheRecourse_)
		{
			if (model->getNumIntegers(1) > 0)
				recourse_[s] = new SolverInterfaceScip(par_);
			else
				recourse_[s] = new DSP_SI(par_);
			recourse_[s]->loadProblem(mat_reco, clbd_reco, cubd_reco, obj_reco, ctype_reco, rlbd_reco, rubd_reco);
			recourse_[s]->setPrintLevel(0);
		}

		/** creating solver interface */
		cglp_[s] = new COIN_OSI;
		cglp_[s]->loadProblem(*mat_reco, clbd_reco, cubd_reco, obj_reco, rlbd_reco, rubd_reco);
		for (int j = 0; j < cglp_[s]->getNumCols(); ++j)
		{
			if (ctype_reco[j] != 'C')
				cglp_[s]->setInteger(j);
		}
		cglp_[s]->messageHandler()->setLogLevel(0);

		/** free memory */
		FREE_MEMORY
	}

	END_TRY_CATCH_RTN(FREE_MEMORY,STO_RTN_ERR)

	FREE_MEMORY;

	return STO_RTN_OK;
#undef FREE_MEMORY
}

#if 0
STO_RTN_CODE TssBdSub::loadProblem(
		int nSubs,
		const double * weight,
		const CoinPackedMatrix ** mat_mp,
		const OsiSolverInterface ** cglp,
		int nAuxvars)
{
	BGN_TRY_CATCH

	nSubs_    = nSubs;
	weight_   = weight;
	mat_mp_   = mat_mp;
	cglp_     = cglp;
	nAuxvars_ = nAuxvars;
	if (cglp_[0])
		solver_infinity_ = cglp_[0]->getInfinity();
	status_ = new int [nSubs_];

#ifndef NO_WARMSTART_TEST
	/** initialize warm start information */
	warm_start_ = new CoinWarmStart * [nSubs_];
	for (int s = nSubs_ - 1; s >= 0; --s)
	{
		warm_start_[s] = NULL;
	}
#endif

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
}
#endif

/** core of generating Benders cut */
void TssBdSub::generateCuts(
		int            ncols,
		const double * x,
		OsiCuts *      cs,
		int            where,
		int            which)
{
#define FREE_MEMORY                    \
	FREE_2D_ARRAY_PTR(nSubs_,Tx);      \
	FREE_2D_ARRAY_PTR(nSubs_, cutval); \
	FREE_ARRAY_PTR(cutrhs);

	if (nSubs_ <= 0)
	{
		DSPdebugMessage("nSubs_ %d\n", nSubs_);
		return;
	}

	double start_time_wall = CoinGetTimeOfDay();

	assert(weight_);
	assert(mat_mp_);
	assert(nAuxvars_ > 0);
	assert(cglp_);

	DSPdebugMessage("\n### START: Bender cut generator ###\n");

	double **     Tx = NULL;
	double ** cutval = NULL; /** dense cut coefficients for each subproblem */
	double *  cutrhs = NULL; /** cut rhs for each subproblem */

	BGN_TRY_CATCH

	/** retrieve the number of rows in subproblem */
	int nrows = 0;
	for (int s = nSubs_ - 1; s >= 0; --s)
	{
		if (!mat_mp_[s]) continue;
		if (excludedScenarios_[s]) continue;
		nrows = mat_mp_[s]->getNumRows();
		break;
	}

#ifdef DSP_DEBUG_DETAIL
	PRINT_ARRAY_MSG(ncols, x, "x");
#endif

	/** allocate memory */
	Tx = new double * [nSubs_];
	for (int s = nSubs_ - 1; s >= 0; --s)
	{
		Tx[s] = NULL;
		if (!mat_mp_[s]) continue;
		if (excludedScenarios_[s]) continue;

		Tx[s] = new double [nrows];
		/** calculate Tx */
		mat_mp_[s]->times(x, Tx[s]);
	}
	cutval = new double * [nSubs_];
	for (int s = nSubs_ - 1; s >= 0; --s)
	{
		cutval[s] = NULL;
		if (!mat_mp_[s]) continue;
		if (excludedScenarios_[s]) continue;

		cutval[s] = new double [ncols];
		CoinZeroN(cutval[s], ncols);
	}
	cutrhs = new double [nSubs_];
	CoinZeroN(cutrhs, nSubs_);

	if (where == TssBd)
	{
		/** loop over subproblems */
		bool doContinue = true;
#ifdef USE_OMP
		/** set number of cores to use */
		omp_set_num_threads(par_->numCores_);
#pragma omp parallel for schedule(dynamic)
#endif
		for (int s = nSubs_ - 1; s >= 0; --s)
		{
			if (!mat_mp_[s]) continue;
			if (excludedScenarios_[s]) continue;
			if (doContinue == false) continue;
			if (which == BothCuts)
				solveOneSubproblem(this, s, x, Tx, cutval, cutrhs);
			else if (which == FeasCut)
				solveOneFeasSubproblem(this, s, x, Tx, cutval, cutrhs);
			else if (which == OptCut)
				solveOneOptSubproblem(this, s, x, Tx, cutval, cutrhs);
			if (status_[s] == STO_STAT_PRIM_INFEASIBLE)
				doContinue = false;
		}
	}
	else if (where == TssDd)
	{
		for (int s = nSubs_ - 1; s >= 0; --s)
		{
			if (!mat_mp_[s]) continue;
			if (excludedScenarios_[s]) continue;
			if (which == BothCuts)
				solveOneSubproblem(this, s, x, Tx, cutval, cutrhs, par_->TssDdAddOptCuts_);
			else if (which == FeasCut)
				solveOneFeasSubproblem(this, s, x, Tx, cutval, cutrhs);
			else if (which == OptCut)
				solveOneOptSubproblem(this, s, x, Tx, cutval, cutrhs);
			if (status_[s] == STO_STAT_PRIM_INFEASIBLE)
				break;
		}
	}

	int status = STO_STAT_OPTIMAL;
	for (int s = nSubs_ - 1; s >= 0; --s)
	{
		if (!mat_mp_[s]) continue;
		if (excludedScenarios_[s]) continue;
		if (status_[s] == STO_STAT_PRIM_INFEASIBLE)
		{
			status = status_[s];
			break;
		}
		else if (status_[s] == STO_STAT_DUAL_INFEASIBLE)
			status = status_[s];
		if (status_[s] == STO_STAT_ABORT ||
			status_[s] == STO_STAT_UNKNOWN)
		{
			status = status_[s];
			break;
		}
	}

	/** construct cuts */
	if (status == STO_STAT_OPTIMAL ||
		 status == STO_STAT_PRIM_INFEASIBLE)
	{
		constructCuts(ncols, cutval, cutrhs, x, cs);
	}

	END_TRY_CATCH(FREE_MEMORY)

	/** collect cut generation time */
	cut_generation_wall_time_ += CoinGetTimeOfDay() - start_time_wall;

	/** free memory */
	FREE_MEMORY

	DSPdebugMessage("### END:   Bender cut generator ###\n\n");
#undef FREE_MEMORY
}

/** solve single scenario recourse */
void TssBdSub::solveSingleRecourse(
		int            scenario, /**< scenario index */
		const double * x,        /**< first-stage variables */
		double &       objval,   /**< second-stage objective value */
		double *       solution  /**< second-stage solution */)
{
#define FREE_MEMORY \
	FREE_ARRAY_PTR(Tx);      \
	FREE_ARRAY_PTR(orgrlbd); \
	FREE_ARRAY_PTR(orgrubd);

	assert(x);
	assert(objval);
	if (par_->TssDdCacheRecourse_)
		assert(recourse_);

	if (!mat_mp_[scenario] || excludedScenarios_[scenario])
		return;

	double * Tx = NULL;
	double * orgrlbd = NULL;
	double * orgrubd = NULL;

	/** retrieve the number of rows in subproblem */
	int nrows = mat_mp_[scenario]->getNumRows();

	/** allocate memory */
	Tx = new double [nrows];
	if (par_->TssDdCacheRecourse_)
	{
		orgrlbd = new double [nrows];
		orgrubd = new double [nrows];
	}

	/** calculate Tx */
	mat_mp_[scenario]->times(x, Tx);

	SolverInterface * si = NULL;
	if (par_->TssDdCacheRecourse_)
	{
		assert(recourse_[scenario]);
		si = recourse_[scenario];
	}
	else
	{
		if (cglp_[scenario]->getNumIntegers() > 0)
		{
			si = new SolverInterfaceScip(par_);
			si->loadProblem(cglp_[scenario], "subprob");
		}
		else
		{
			si = new DSP_SI(par_, cglp_[scenario]);
		}
		si->setPrintLevel(cglp_[scenario]->messageHandler()->logLevel());
	}

	/** local variables */
	const double * rlbd = si->getRowLower();
	const double * rubd = si->getRowUpper();

	/** copy original row bounds */
	if (par_->TssDdCacheRecourse_)
	{
		CoinCopyN(rlbd, si->getNumRows(), orgrlbd);
		CoinCopyN(rubd, si->getNumRows(), orgrubd);
	}

	/** update row bounds */
	for (int i = nrows - 1; i >= 0; --i)
	{
		if (rlbd[i] > -COIN_DBL_MAX)
			si->setRowLower(i, rlbd[i] - Tx[i]);
		if (rubd[i] < COIN_DBL_MAX)
			si->setRowUpper(i, rubd[i] - Tx[i]);
	}

	/** solve */
	si->solve();

	/** solution status */
	status_[scenario] = si->getStatus();

	/** collect solution */
	switch (status_[scenario])
	{
	case STO_STAT_OPTIMAL:
	case STO_STAT_STOPPED_GAP:
	{
		/** objective */
		objval = si->getPrimalBound();
		/** solution */
		if (solution)
			CoinCopyN(si->getSolution(), si->getNumCols(), solution);
		break;
	}
	default:
		objval = COIN_DBL_MAX;
		break;
	}

	/** restore row bounds */
	if (par_->TssDdCacheRecourse_)
	{
		for (int i = nrows - 1; i >= 0; --i)
		{
			si->setRowLower(i, orgrlbd[i]);
			si->setRowUpper(i, orgrubd[i]);
		}
		si = NULL;
	}
	else
		FREE_PTR(si);

	/** null */
	rlbd = NULL;
	rubd = NULL;

	/** free memory */
	FREE_MEMORY;

#undef FREE_MEMORY
}

/** solve recourse (without cut generation) */
void TssBdSub::solveRecourse(
		const double * x,        /**< first-stage variables */
		double *       objval,   /**< second-stage objective values */
		double **      solution, /**< second-stage solutions */
		int            ncores    /**< number of cores used to run in parallel */)
{
#define FREE_MEMORY \
	FREE_2D_ARRAY_PTR(ncores,Tx);      \
	FREE_2D_ARRAY_PTR(ncores,orgrlbd); \
	FREE_2D_ARRAY_PTR(ncores,orgrubd);

	assert(x);
	assert(objval);
	if (par_->TssDdCacheRecourse_)
		assert(recourse_);

	double ** Tx = NULL;
	double ** orgrlbd = NULL;
	double ** orgrubd = NULL;

	/** retrieve the number of rows in subproblem */
	int nrows = 0;
	for (int s = nSubs_ - 1; s >= 0; --s)
	{
		if (!mat_mp_[s] || excludedScenarios_[s]) continue;
		nrows = mat_mp_[s]->getNumRows();
		break;
	}

	/** allocate memory */
	Tx = new double * [ncores];
	for (int s = ncores - 1; s >= 0; --s)
	{
		Tx[s] = NULL;
		if (!mat_mp_[s] || excludedScenarios_[s]) continue;

		Tx[s] = new double [nrows];
	}
	if (par_->TssDdCacheRecourse_)
	{
		orgrlbd = new double * [ncores];
		orgrubd = new double * [ncores];
		for (int s = ncores - 1; s >= 0; --s)
		{
			orgrlbd[s] = new double [nrows];
			orgrubd[s] = new double [nrows];
		}
	}

#ifdef USE_OMP
	/** set number of cores to use */
	omp_set_num_threads(ncores);
#pragma omp parallel for schedule(dynamic) default(shared)
#endif
	for (int s = nSubs_ - 1; s >= 0; --s)
	{
		if (!mat_mp_[s] || excludedScenarios_[s]) continue;

#ifdef USE_OMP
		/** get thread ID */
		int tid = omp_get_thread_num();
#else
		int tid = 0;
#endif

		/** calculate Tx */
		mat_mp_[s]->times(x, Tx[tid]);

		SolverInterface * si = NULL;
		if (par_->TssDdCacheRecourse_)
		{
			assert(recourse_[s]);
			si = recourse_[s];
		}
		else
		{
			if (cglp_[s]->getNumIntegers() > 0)
			{
				si = new SolverInterfaceScip(par_);
				si->loadProblem(cglp_[s], "subprob");
			}
			else
			{
				si = new DSP_SI(par_, cglp_[s]);
			}
			si->setPrintLevel(cglp_[s]->messageHandler()->logLevel());
		}

		/** local variables */
		const double * rlbd = si->getRowLower();
		const double * rubd = si->getRowUpper();

		/** copy original row bounds */
		if (par_->TssDdCacheRecourse_)
		{
			CoinCopyN(rlbd, si->getNumRows(), orgrlbd[tid]);
			CoinCopyN(rubd, si->getNumRows(), orgrubd[tid]);
		}

		/** update row bounds */
		for (int i = nrows - 1; i >= 0; --i)
		{
			if (rlbd[i] > -COIN_DBL_MAX)
				si->setRowLower(i, rlbd[i] - Tx[tid][i]);
			if (rubd[i] < COIN_DBL_MAX)
				si->setRowUpper(i, rubd[i] - Tx[tid][i]);
			//printf("row %d: rlbd %f rubd %f Tx %f\n", i, rlbd[i], rubd[i], Tx[tid][i]);
		}

#if 0
		char filename[128];
		sprintf(filename, "ub%d.mps", s);
		si->writeMps(filename);
#endif
		/** solve */
		si->solve();

		/** solution status */
		status_[s] = si->getStatus();

		/** collect solution */
		switch (status_[s])
		{
		case STO_STAT_OPTIMAL:
		case STO_STAT_STOPPED_GAP:
		{
			/** objective */
			objval[s] = si->getPrimalBound();
			/** solution */
			if (solution)
				CoinCopyN(si->getSolution(), si->getNumCols(), solution[s]);
			break;
		}
		default:
			objval[s] = COIN_DBL_MAX;
			break;
		}

		/** restore row bounds */
		if (par_->TssDdCacheRecourse_)
		{
			for (int i = nrows - 1; i >= 0; --i)
			{
				si->setRowLower(i, orgrlbd[tid][i]);
				si->setRowUpper(i, orgrubd[tid][i]);
			}
			si = NULL;
		}
		else
			FREE_PTR(si);

		/** null */
		rlbd = NULL;
		rubd = NULL;
	}

	/** free memory */
	FREE_MEMORY;

#undef FREE_MEMORY
}

/** solve one subproblem. this is a body of loop in gutsOfGenerateCuts */
void TssBdSub::solveOneSubproblem(
		TssBdSub *     cgl,
		int            s,            /**< scenario index */
		const double * x,            /**< first-stage solution */
		double **      Tx,           /**< Tx */
		double **      cutval,       /**< Benders cut body */
		double *       cutrhs,       /**< Benders cut RHS */
		int            enableOptCuts /**< whether to generate optimality cuts or not */)
{
	DSPdebugMessage("Scenario %d\n", s);

	/** local variables */
	const double * rlbd = cgl->cglp_[s]->getRowLower();
	const double * rubd = cgl->cglp_[s]->getRowUpper();
	const double * clbd = cgl->cglp_[s]->getColLower();
	const double * cubd = cgl->cglp_[s]->getColUpper();
	const double * pi   = NULL; /** dual variables */
	const double * rc   = NULL; /** reduced costs */

	/** clone CGLP */
#if 0
	OsiSolverInterface * cglp = cgl->cglp_[s]->clone();
#else
	OsiSolverInterface * cglp = new COIN_OSI;
	cglp->loadProblem(*(cgl->cglp_[s]->getMatrixByCol()),
			cgl->cglp_[s]->getColLower(), cgl->cglp_[s]->getColUpper(),
			cgl->cglp_[s]->getObjCoefficients(),
			cgl->cglp_[s]->getRowLower(), cgl->cglp_[s]->getRowUpper());
	cglp->setWarmStart(cgl->warm_start_[s]);
#endif

	int nrows = cglp->getNumRows();

	/** loop over CGLP rows to update row bounds */
	for (int i = nrows - 1; i >= 0; --i)
	{
		if (rlbd[i] > -COIN_DBL_MAX)
			cglp->setRowLower(i, rlbd[i] - Tx[s][i]);
		if (rubd[i] < COIN_DBL_MAX)
			cglp->setRowUpper(i, rubd[i] - Tx[s][i]);
	}

#ifdef SOLVER_FEASIBILITY_FIRST
	/** solve feasibility problem */
#if 0
#if 0
	OsiSolverInterface * si_copy = cglp->clone();
#else
	OsiSolverInterface * si_copy = new COIN_OSI;
	si_copy->loadProblem(*(cglp->getMatrixByCol()),
			cglp->getColLower(), cglp->getColUpper(),
			cglp->getObjCoefficients(),
			cglp->getRowLower(), cglp->getRowUpper());
	si_copy->setWarmStart(cglp->getWarmStart());
#endif
	int nAddedCols = 0;
	solveFeasProblem(si_copy, nAddedCols);

	/** get warmstart */
	if (!cgl->warm_start_[s])
		cgl->warm_start_[s] = si_copy->getWarmStart();

	/** infeasible? */
	if (si_copy->getObjValue() > 1.e-6)
	{
		/** save solution status */
		cgl->status_[s] = STO_STAT_PRIM_INFEASIBLE;

		/** calculate feasibility cut elements */
		calculateCutElements(
				si_copy->getNumRows(), si_copy->getNumCols(), cgl->mat_mp_[s],
				rlbd, rubd, si_copy->getColLower(), si_copy->getColUpper(),
				si_copy->getRowPrice(), si_copy->getReducedCost(), cutval[s], cutrhs[s]);

		/** free memory */
		FREE_PTR(si_copy);
		FREE_PTR(cglp);

		return;
	}
	FREE_PTR(si_copy);

	if (!enableOptCuts)
	{
		/** free memory */
		FREE_PTR(cglp);
		return;
	}
#else
	int nAddedCols = 0;

	/** change to and solve feasibility problem */
	solveFeasProblem(cglp, nAddedCols);

	/** get warmstart */
	if (!cgl->warm_start_[s])
	{
		cgl->warm_start_[s] = cglp->getWarmStart();
		CoinWarmStartBasis * basis = dynamic_cast<CoinWarmStartBasis*>(cgl->warm_start_[s]);
		basis->resize(cglp->getNumRows(), cglp->getNumCols() - nAddedCols);
		basis = NULL;
	}

	/** infeasible? */
	if (cglp->getObjValue() > 1.e-6)
	{
		/** save solution status */
		cgl->status_[s] = STO_STAT_PRIM_INFEASIBLE;

		/** calculate feasibility cut elements */
		calculateCutElements(
				cglp->getNumRows(), cglp->getNumCols(), cgl->mat_mp_[s],
				rlbd, rubd, cglp->getColLower(), cglp->getColUpper(),
				cglp->getRowPrice(), cglp->getReducedCost(), cutval[s], cutrhs[s]);

		/** free memory */
		FREE_PTR(cglp);

		return;
	}

	if (!enableOptCuts)
	{
		cgl->status_[s] = STO_STAT_FEASIBLE;
		/** free memory */
		FREE_PTR(cglp);
		return;
	}
	/** change to original problem */
	chgToOrgProblem(cglp, cgl->cglp_[s]->getObjCoefficients(), nAddedCols);
#endif
#endif

	//printf("<<< Optimality Problem >>>\n");
	cglp->messageHandler()->setLogLevel(0);
	cglp->setHintParam(OsiDoScale);
	cglp->setHintParam(OsiDoPresolveInResolve);
	cglp->setHintParam(OsiDoDualInResolve);

	/** set warmstart */
	cglp->setWarmStart(cgl->warm_start_[s]);

	/** solve */
	cglp->resolve();

	/** update warmstart */
	FREE_PTR(cgl->warm_start_[s]);
	cgl->warm_start_[s] = cglp->getWarmStart();

	DSPdebugMessage("  objective value %E\n", cglp->getObjValue());

	/** solution status */
	if (cglp->isProvenOptimal() ||
		cglp->isIterationLimitReached())
	{
		/** save solution status */
		cgl->status_[s] = STO_STAT_OPTIMAL;
		DSPdebugMessage("  solution status: optimal\n");

		/** TODO: add parametric cuts */

		/** get dual variables and reduced costs */
		pi = cglp->getRowPrice();
		rc = cglp->getReducedCost();
#ifdef DSP_DEBUG2
		for (int i = 0; i < cglp->getNumRows(); ++i)
		{
			printf("  rlbd[%d] = %E, rubd[%d] = %E, pi[%d] = %E\n",
					i, rlbd[i], i, rubd[i], i, pi[i]);
		}
		for (int j = 0; j < cglp->getNumCols(); ++j)
		{
			printf("  y[%d] = %E, rc[%d] = %E\n", j, cglp->getColSolution()[j], j, rc[j]);
		}
#endif

		/** calculate cut elements */
		calculateCutElements(cglp->getNumRows(), cglp->getNumCols(),
				cgl->mat_mp_[s], rlbd, rubd, clbd, cubd, pi, rc, cutval[s], cutrhs[s]);
	}
	else if (cglp->isProvenPrimalInfeasible())
	{
		/** save solution status */
		cgl->status_[s] = STO_STAT_PRIM_INFEASIBLE;
		DSPdebugMessage("  solution status: primal infeasible\n");
#ifndef SOLVER_FEASIBILITY_FIRST
		/** solve feasibility problem */
#if 0
		OsiSolverInterface * si_copy = cglp->clone();
#else
		OsiSolverInterface * si_copy = new COIN_OSI;
		si_copy->loadProblem(*(cglp->getMatrixByCol()),
				cglp->getColLower(), cglp->getColUpper(),
				cglp->getObjCoefficients(),
				cglp->getRowLower(), cglp->getRowUpper());
		si_copy->setWarmStart(cglp->getWarmStart());
#endif
		solveFeasProblem(si_copy);

		/** calculate feasibility cut elements */
		calculateCutElements(
				si_copy->getNumRows(), si_copy->getNumCols(), cgl->mat_mp_[s],
				rlbd, rubd, si_copy->getColLower(), si_copy->getColUpper(),
				si_copy->getRowPrice(), si_copy->getReducedCost(), cutval[s], cutrhs[s]);

		/** free memory */
		FREE_PTR(si_copy);
#endif
	}
	else if (cglp->isProvenDualInfeasible())
	{
		/** save solution status */
		cgl->status_[s] = STO_STAT_DUAL_INFEASIBLE;
		DSPdebugMessage("  solution status: dual infeasible\n");
	}
	else if (cglp->isAbandoned() ||
			 cglp->isPrimalObjectiveLimitReached() ||
			 cglp->isDualObjectiveLimitReached())
	{
		cgl->status_[s] = STO_STAT_STOPPED_UNKNOWN;
		DSPdebugMessage("  solution status: stopped unknown\n");
	}
	else
	{
		cgl->status_[s] = STO_STAT_UNKNOWN;
		DSPdebugMessage("  solution status: unknown\n");
	}

	/** free clone */
	FREE_PTR(cglp);
}

/** solve one feasibility subproblem. this is a body of loop in gutsOfGenerateCuts */
void TssBdSub::solveOneFeasSubproblem(
		TssBdSub *     cgl,
		int            s,      /**< scenario index */
		const double * x,      /**< first-stage solution */
		double **      Tx,     /**< Tx */
		double **      cutval, /**< Benders cut body */
		double *       cutrhs  /**< Benders cut RHS */)
{
	DSPdebugMessage("Scenario %d\n", s);

	/** local variables */
	const double * rlbd = cgl->cglp_[s]->getRowLower();
	const double * rubd = cgl->cglp_[s]->getRowUpper();

	/** clone CGLP */
#if 0
	OsiSolverInterface * cglp = cgl->cglp_[s]->clone();
#else
	OsiSolverInterface * cglp = new COIN_OSI;
	cglp->loadProblem(*(cgl->cglp_[s]->getMatrixByCol()),
			cgl->cglp_[s]->getColLower(), cgl->cglp_[s]->getColUpper(),
			cgl->cglp_[s]->getObjCoefficients(),
			cgl->cglp_[s]->getRowLower(), cgl->cglp_[s]->getRowUpper());
	cglp->setWarmStart(cgl->warm_start_[s]);
#endif

	int nrows = cglp->getNumRows();

	/** loop over CGLP rows to update row bounds */
	for (int i = nrows - 1; i >= 0; --i)
	{
		if (rlbd[i] > -COIN_DBL_MAX)
			cglp->setRowLower(i, rlbd[i] - Tx[s][i]);
		if (rubd[i] < COIN_DBL_MAX)
			cglp->setRowUpper(i, rubd[i] - Tx[s][i]);
	}

	/** solve feasibility problem */
	int nAddedCols = 0;

	/** change to and solve feasibility problem */
	solveFeasProblem(cglp, nAddedCols);

	/** get warmstart */
	if (!cgl->warm_start_[s])
	{
		cgl->warm_start_[s] = cglp->getWarmStart();
		CoinWarmStartBasis * basis = dynamic_cast<CoinWarmStartBasis*>(cgl->warm_start_[s]);
		basis->resize(cglp->getNumRows(), cglp->getNumCols() - nAddedCols);
		basis = NULL;
	}

	/** infeasible? */
	if (cglp->getObjValue() > 1.e-6)
	{
		/** save solution status */
		cgl->status_[s] = STO_STAT_PRIM_INFEASIBLE;

		/** calculate feasibility cut elements */
		calculateCutElements(
				cglp->getNumRows(), cglp->getNumCols(), cgl->mat_mp_[s],
				rlbd, rubd, cglp->getColLower(), cglp->getColUpper(),
				cglp->getRowPrice(), cglp->getReducedCost(), cutval[s], cutrhs[s]);

		/** free memory */
		FREE_PTR(cglp);

		return;
	}

	cgl->status_[s] = STO_STAT_FEASIBLE;

	/** free memory */
	FREE_PTR(cglp);

	return;
}

/** solve one feasibility subproblem. this is a body of loop in gutsOfGenerateCuts */
void TssBdSub::solveOneOptSubproblem(
		TssBdSub *     cgl,
		int            s,      /**< scenario index */
		const double * x,      /**< first-stage solution */
		double **      Tx,     /**< Tx */
		double **      cutval, /**< Benders cut body */
		double *       cutrhs  /**< Benders cut RHS */)
{
	DSPdebugMessage("Scenario %d\n", s);

	/** local variables */
	const double * rlbd = cgl->cglp_[s]->getRowLower();
	const double * rubd = cgl->cglp_[s]->getRowUpper();
	const double * clbd = cgl->cglp_[s]->getColLower();
	const double * cubd = cgl->cglp_[s]->getColUpper();
	const double * pi   = NULL; /** dual variables */
	const double * rc   = NULL; /** reduced costs */

	/** clone CGLP */
#if 0
	OsiSolverInterface * cglp = cgl->cglp_[s]->clone();
#else
	OsiSolverInterface * cglp = new COIN_OSI;
	cglp->loadProblem(*(cgl->cglp_[s]->getMatrixByCol()),
			cgl->cglp_[s]->getColLower(), cgl->cglp_[s]->getColUpper(),
			cgl->cglp_[s]->getObjCoefficients(),
			cgl->cglp_[s]->getRowLower(), cgl->cglp_[s]->getRowUpper());
	cglp->setWarmStart(cgl->warm_start_[s]);
#endif

	int nrows = cglp->getNumRows();

	/** loop over CGLP rows to update row bounds */
	for (int i = nrows - 1; i >= 0; --i)
	{
		if (rlbd[i] > -COIN_DBL_MAX)
			cglp->setRowLower(i, rlbd[i] - Tx[s][i]);
		if (rubd[i] < COIN_DBL_MAX)
			cglp->setRowUpper(i, rubd[i] - Tx[s][i]);
	}

	//printf("<<< Optimality Problem >>>\n");
	cglp->messageHandler()->setLogLevel(0);
	cglp->setHintParam(OsiDoScale);
	cglp->setHintParam(OsiDoPresolveInResolve);
	cglp->setHintParam(OsiDoDualInResolve);

	/** set warmstart */
	cglp->setWarmStart(cgl->warm_start_[s]);

	/** solve */
	cglp->resolve();

	/** update warmstart */
	FREE_PTR(cgl->warm_start_[s]);
	cgl->warm_start_[s] = cglp->getWarmStart();

	DSPdebugMessage("  objective value %E\n", cglp->getObjValue());

	/** solution status */
	if (cglp->isProvenOptimal() ||
		cglp->isIterationLimitReached())
	{
		/** save solution status */
		cgl->status_[s] = STO_STAT_OPTIMAL;
		DSPdebugMessage("  solution status: optimal\n");

		/** TODO: add parametric cuts */

		/** get dual variables and reduced costs */
		pi = cglp->getRowPrice();
		rc = cglp->getReducedCost();

		/** calculate cut elements */
		calculateCutElements(cglp->getNumRows(), cglp->getNumCols(),
				cgl->mat_mp_[s], rlbd, rubd, clbd, cubd, pi, rc, cutval[s], cutrhs[s]);
	}
	else if (cglp->isProvenPrimalInfeasible())
	{
		/** save solution status */
		cgl->status_[s] = STO_STAT_PRIM_INFEASIBLE;
		DSPdebugMessage("  solution status: primal infeasible\n");
	}
	else if (cglp->isProvenDualInfeasible())
	{
		/** save solution status */
		cgl->status_[s] = STO_STAT_DUAL_INFEASIBLE;
		DSPdebugMessage("  solution status: dual infeasible\n");
	}
	else if (cglp->isAbandoned() ||
			 cglp->isPrimalObjectiveLimitReached() ||
			 cglp->isDualObjectiveLimitReached())
	{
		cgl->status_[s] = STO_STAT_STOPPED_UNKNOWN;
		DSPdebugMessage("  solution status: stopped unknown\n");
	}
	else
	{
		cgl->status_[s] = STO_STAT_UNKNOWN;
		DSPdebugMessage("  solution status: unknown\n");
	}

	/** free clone */
	FREE_PTR(cglp);

	return;
}

/** solve feasibility problem */
int TssBdSub::solveFeasProblem(
		OsiSolverInterface * si, /**< [in] subproblem solver interface */
		int & nAddedCols         /**< [out] number of columns added */)
{
	BGN_TRY_CATCH

	for (int j = 0; j < si->getNumCols(); ++j)
		si->setObjCoeff(j, 0.0);
	std::vector<int> columnStarts, rows;
	std::vector<double> elements, clbd, cubd, obj;
	for (int i = 0; i < si->getNumRows(); ++i)
	{
		columnStarts.push_back(elements.size());
		rows.push_back(i);
		if (si->getRowSense()[i] == 'G')
			elements.push_back(1.0);
		else
			elements.push_back(-1.0);
		clbd.push_back(0.0);
		cubd.push_back(COIN_DBL_MAX);
		obj.push_back(1.0);
		if (si->getRowSense()[i] == 'E' ||
			si->getRowSense()[i] == 'R')
		{
			columnStarts.push_back(elements.size());
			rows.push_back(i);
			elements.push_back(1.0);
			clbd.push_back(0.0);
			cubd.push_back(COIN_DBL_MAX);
			obj.push_back(1.0);
		}
	}
	columnStarts.push_back(elements.size());
	nAddedCols = columnStarts.size() - 1;
	si->addCols(nAddedCols, &columnStarts[0], &rows[0], &elements[0], &clbd[0], &cubd[0], &obj[0]);

#ifdef DSP_DEBUG
	si->messageHandler()->setLogLevel(1);
#else
	si->messageHandler()->setLogLevel(0);
#endif
	si->setHintParam(OsiDoScale);
	si->setHintParam(OsiDoPresolveInInitial);
	si->setHintParam(OsiDoDualInInitial);

#if 0
	if (si->getWarmStart())
	{
		/** set basis */
		CoinWarmStartBasis * basis = dynamic_cast<CoinWarmStartBasis*>(si->getWarmStart());
		basis->resize(si->getNumRows(), si->getNumCols());
		for (int i = 0; i < nAddedCols; ++i)
			basis->setStructStatus(si->getNumCols() - nAddedCols + i, CoinWarmStartBasis::basic);
		si->setWarmStart(basis);

		si->resolve();

		FREE_PTR(basis);
	}
	else
#endif
	{
		si->initialSolve();
	}
	//si_copy->writeMps("debug");

	/** should be always optimal */
	assert(si->isProvenOptimal());

	END_TRY_CATCH_RTN(;, 1)

	return 0;
}

/** change feasibility problem to optimality problem */
int TssBdSub::chgToOrgProblem(
		OsiSolverInterface * si, /**< [in] subproblem solver interface */
		const double * obj,      /**< [in] original objective function */
		int & nAddedCols         /**< [out] number of columns added */)
{
#define FREE_MEMORY FREE_ARRAY_PTR(delind)

	int * delind = NULL;

	BGN_TRY_CATCH

	/** delete columns */
	delind = new int [nAddedCols];
	CoinIotaN(delind, nAddedCols, si->getNumCols() - nAddedCols);
	si->deleteCols(nAddedCols, delind);

	/** set original objective */
	si->setObjective(obj);

	END_TRY_CATCH_RTN(FREE_MEMORY;, 1)

	FREE_MEMORY;

	return 0;
#undef FREE_MEMORY
}

/** construct cuts (e.g. aggregation) */
/** TODO: Only supports a simple aggregation for now. */
int TssBdSub::constructCuts(
		int            ncols,         /**< [in] number of columns in the master */
		double **      values,        /**< [in] dense vector of cut coefficients in the size of nSubs_ */
		double *       rhs,           /**< [in] array of cut rhs in the size of nSubs_ */
		const double * x,             /**< [in] current solution to calculate effectiveness */
		OsiCuts *      cuts           /**< [out] a set of cuts generated */)
{
#define FREE_MEMORY                       \
	FREE_2D_ARRAY_PTR(nAuxvars_, aggval); \
	FREE_ARRAY_PTR(aggrhs);

	int s, j;
	double ** aggval  = NULL; /** aggregated dense cut coefficients */
	double *  aggrhs  = NULL; /** aggregated cut rhs */
	CoinPackedVector vec;

	BGN_TRY_CATCH

	/** allocate memory */
	aggval = new double * [nAuxvars_];
	for (s = nAuxvars_ - 1; s >= 0; --s)
	{
		aggval[s] = new double [ncols];
		CoinZeroN(aggval[s], ncols);
	}
	aggrhs = new double [nAuxvars_];
	CoinZeroN(aggrhs, nAuxvars_);

	for (s = nSubs_ - 1; s >= 0; --s)
	{
		if (!mat_mp_[s]) continue;
		if (excludedScenarios_[s]) continue;
		/** feasibility cut? */
		if (status_[s] == STO_STAT_PRIM_INFEASIBLE)
		{
			/** initialize vector */
			vec.clear();

			/** set it as sparse */
			for (j = 0; j < ncols; ++j)
			{
				if (fabs(values[s][j]) > 1E-10)
				{
					vec.insert(j, weight_[s] * values[s][j]);
				}
			}

			rhs[s] *= weight_[s];
			if (fabs(rhs[s]) < 1E-10)
				rhs[s] = 0.0;

			OsiRowCut rc;
			rc.setRow(vec);
			rc.setUb(COIN_DBL_MAX); /** for minimization */
			rc.setLb(rhs[s]);

			DSPdebug(rc.print());
			cuts->insert(rc);

			break;
		}

		/** consider only optimality cuts */
		if (status_[s] != STO_STAT_OPTIMAL) continue;

		/** calculate weighted aggregation of cuts */
		int ind_aux = s % nAuxvars_;
		for (j = ncols - 1; j >= 0; --j)
		{
			aggval[ind_aux][j] += weight_[s] * values[s][j];
			DSPdebugMessage("  aggval[%d][%d] %E weight_[%d] %E cutval[%d][%d] %E\n",
					ind_aux, j, aggval[ind_aux][j],
					s, weight_[s],
					s, j, values[s][j]);
		}
		aggrhs[ind_aux] += weight_[s] * rhs[s];
		DSPdebugMessage("  aggval[%d][%d] %E weight_[%d] %E\n",
				ind_aux, ncols - nAuxvars_ + ind_aux, aggval[ind_aux][ncols - nAuxvars_ + ind_aux],
				s, weight_[s]);
		DSPdebugMessage("  aggrhs[%d] %E weight_[%d] %E cutrhs[%d] %E\n",
				ind_aux, aggrhs[ind_aux],
				s, weight_[s],
				s, rhs[s]);
	}

	/** construct cuts to pass */
	for (s = nAuxvars_ - 1; s >= 0; --s)
	{
		/** auxiliary variable coefficient */
		aggval[s][ncols - nAuxvars_ + s] = 1;

		/** initialize vector */
		vec.clear();

		/** set it as sparse */
		for (j = 0; j < ncols; ++j)
		{
			if (fabs(aggval[s][j]) > 1E-10)
			{
				vec.insert(j, aggval[s][j]);
			}
		}

		if (fabs(aggrhs[s]) < 1E-10)
			aggrhs[s] = 0.0;

		/** effective? */
		//if (vec.getNumElements() > 1)
		{
			OsiRowCut rc;
			rc.setRow(vec);
			rc.setUb(COIN_DBL_MAX); /** for minimization */
			rc.setLb(aggrhs[s]);

			DSPdebug(rc.print());
			cuts->insert(rc);
		}
	}

	END_TRY_CATCH_RTN(FREE_MEMORY,1)

	FREE_MEMORY;

	return 0;
#undef FREE_MEMORY
}

/** calculate cut elements and rhs */
int TssBdSub::calculateCutElements(
		int nrows,                         /**< [in] number of rows in subproblem */
		int ncols,                         /**< [in] number of columns in subproblem */
		const CoinPackedMatrix * mat_tech, /**< [in] technology matrix */
		const double * rlbd,               /**< [in] row lower bounds */
		const double * rubd,               /**< [in] row upper bounds */
		const double * clbd,               /**< [in] column lower bounds */
		const double * cubd,               /**< [in] column upper bounds */
		const double * pi,                 /**< [in] dual variables corresponding to constraints */
		const double * rc,                 /**< [in] reduced cost corresponding to variables */
		double *       cutval,             /**< [out] cut coefficients */
		double &       cutrhs              /**< [out] cut rhs */)
{
	int i, j;

	/** calculate pi^T T */
	mat_tech->transposeTimes(pi, cutval);

	/** calculate pi^T h */
	/** loop over rows */
	for (i = nrows - 1; i >= 0; --i)
	{
		//if (fabs(pi[i]) < 1.0e-10) continue;
        if (rubd[i] >= COIN_DBL_MAX)
        {
            /** Ax >= b */
        	cutrhs += pi[i] * rlbd[i];
        }
        else if (rlbd[i] <= -COIN_DBL_MAX)
        {
            /** Ax <= b */
        	cutrhs += pi[i] * rubd[i];
        }
        else if (rlbd[i] == rubd[i])
        {
            /** Ax = b */
        	cutrhs += pi[i] * rlbd[i];
        }
        else if (pi[i] >= 0.0)
        {
            /** l <= Ax <= u, bounded by l (-y + s_L - s_U = 0) */
        	cutrhs += pi[i] * rlbd[i];
        }
        else
        {
            /** l <= Ax <= u, bounded by u */
        	cutrhs += pi[i] * rubd[i];
        }
		//printf("cutrhs %e, pi %e rlbd %e rubd %e\n", cutrhs, pi[i], rlbd[i], rubd[i]);
	}
//	printf("cutrhs %E\n", cutrhs);

	/** loop over columns */
	for (j = ncols - 1; j >= 0; --j)
	{
		//if (fabs(rc[j]) < 1.0e-10) continue;
		//cutval[j] += rc[j];
        if (cubd[j] >= COIN_DBL_MAX)
        {
            /** x_j >= l_j */
            cutrhs += rc[j] * clbd[j];
        }
        else if (clbd[j] <= -COIN_DBL_MAX)
        {
            /** x_j <= u_j */
            cutrhs += rc[j] * cubd[j];
        }
        else if (clbd[j] == cubd[j])
        {
            /** l_j = x_j = u_j */
            cutrhs += rc[j] * clbd[j];
        }
        else if (rc[j] >= 0.0)
        {
            /** l_j <= x_j <= u_j, bounded by l */
            cutrhs += rc[j] * clbd[j];
        }
        else
        {
            /** l_j <= x_j <= u_j, bounded by u */
            cutrhs += rc[j] * cubd[j];
        }
		//printf("cutrhs %e, rc %e clbd %e cubd %e\n", cutrhs, rc[j], clbd[j], cubd[j]);
	}

//	printf("cutrhs %E\n", cutrhs);

	return 0;
}

