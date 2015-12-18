/*
 * TssBd.cpp
 *
 *  Created on: Dec 5, 2014
 *      Author: kibaekkim
 */

//#define DSP_DEBUG

/** DSP */
#include "Solver/TssBd.h"
#include "SolverInterface/SCIPconshdlrBenders.h"
#include "SolverInterface/SolverInterfaceSpx.h"
#include "SolverInterface/SolverInterfaceClp.h"
#include "SolverInterface/SolverInterfaceScip.h"
#include "Utility/StoMessage.h"

TssBd::TssBd() :
	TssSolver(),
	si_(NULL),
	naugs_(0),
	augs_(NULL),
	naux_(0),
	obj_aux_(NULL),
	clbd_aux_(NULL),
	cubd_aux_(NULL),
	parNumCores_(1),
	parCutPriority_(-200000),
	parInitLbAlgo_(SEPARATE_LP)
{
	/** nothing to do */
}

TssBd::~TssBd()
{
	FREE_PTR(si_);
	FREE_ARRAY_PTR(augs_);
	FREE_ARRAY_PTR(obj_aux_);
	FREE_ARRAY_PTR(clbd_aux_);
	FREE_ARRAY_PTR(cubd_aux_);
	naugs_ = 0;
	naux_ = 0;
}

/** initialize */
STO_RTN_CODE TssBd::initialize()
{
	/** parameters */
	parNumCores_    = par_->getIntParam("BD/NUM_CORES");
	parCutPriority_ = par_->getIntParam("BD/CUT_PRIORITY");
	parInitLbAlgo_  = par_->getIntParam("BD/INIT_LB_ALGO");

	/** set objective coefficients for auxiliary variables */
	CoinZeroN(obj_aux_, naux_);
	for (int s = 0; s < model_->getNumScenarios(); ++s)
		obj_aux_[s % naux_] += model_->getProbability()[s];

	return STO_RTN_OK;
}

/** solve */
STO_RTN_CODE TssBd::solve()
{
#define FREE_MEMORY    \
	FREE_PTR(tssbdsub)

	assert(model_);
	assert(par_);

	double tic;
	TssBdSub * tssbdsub = NULL; /**< cut generator */
	double lowerbound = 0.0;

	BGN_TRY_CATCH

	/** initialize */
	initialize();

	/** solution time */
	double swtime = CoinGetTimeOfDay();

	message_->print(1, "\nPhase 1:\n");
	message_->print(1, "Creating Benders sub problems ...");
	tic = CoinGetTimeOfDay();

	/** configure Benders cut generator */
	/** This does NOT make deep copies. So, DO NOT RELEASE POINTERS. */
	tssbdsub = new TssBdSub(par_);
	for (int s = 0; s < model_->getNumScenarios(); ++s)
		tssbdsub->scenarios_.push_back(s);
	STO_RTN_CHECK_THROW(tssbdsub->loadProblem(model_, naugs_, augs_, naux_), "loadProblem", "TssBdSub");

	message_->print(1, " (%.2f sec)\n", CoinGetTimeOfDay() - tic);

	/** solution time */
	double stime = clockType_ ? CoinGetTimeOfDay() : CoinCpuTime();

	message_->print(1, "Finding global lower bound ...");
	tic = CoinGetTimeOfDay();

	/** find lower bound */
	/** TODO: replace this with TssDd */
	STO_RTN_CHECK_THROW(findLowerBound(lowerbound), "findLowerBound", "TssBd");

	message_->print(1, " (%.2f sec) -> Lower bound %e\n", CoinGetTimeOfDay() - tic, lowerbound);
	message_->print(1, "Creating master problem instance ...");
	tic = CoinGetTimeOfDay();

	/** construct master problem */
	STO_RTN_CHECK_THROW(constructMasterProblem(tssbdsub, lowerbound), "constructMasterProblem", "TssBd");

	message_->print(1, " (%.2f sec)\n", CoinGetTimeOfDay() - tic);
	message_->print(1, "\nPhase 2:\n");

	time_remains_ -= CoinGetTimeOfDay() - swtime;

	/** use L-shaped method to find lower bound */
	if (model_->getNumCoreIntegers() == 0)
	{
		/** configure LP */
		STO_RTN_CHECK_THROW(configureSLP(), "configureSLP", "TssBd");

		/** solve LP */
		STO_RTN_CHECK_THROW(solveSLP(tssbdsub), "solveSLP", "TssBd");
	}
	else
	{
		/** configure MILP */
		STO_RTN_CHECK_THROW(configureSMILP(tssbdsub), "configureSMILP", "TssBd");

		/** solve MILP */
		STO_RTN_CHECK_THROW(solveSMILP(), "solveSMILP", "TssBd");
	}

	/** solution time */
	solutionTime_ = (clockType_ ? CoinGetTimeOfDay() : CoinCpuTime()) - stime;

	message_->print(1, "\nCollecting results ...");
	tic = CoinGetTimeOfDay();

	/** collect solution */
	switch (status_)
	{
	case STO_STAT_OPTIMAL:
	case STO_STAT_LIM_ITERorTIME:
	case STO_STAT_STOPPED_GAP:
	case STO_STAT_STOPPED_NODE:
	case STO_STAT_STOPPED_TIME:
		{
			const double * solution = si_->getSolution();
			if (solution)
			{
				/** first-stage solution */
				assert(solution_);
				CoinCopyN(solution, model_->getNumCols(0), solution_);

				/** second-stage solution */
				double * objval_reco = new double [model_->getNumScenarios()];
				double ** solution_reco = new double * [model_->getNumScenarios()];
				for (int s = 0; s < model_->getNumScenarios(); ++s)
					solution_reco[s] = new double [model_->getNumCols(1)];

				for (int s = 0; s < naugs_; ++s)
					CoinCopyN(solution + model_->getNumCols(0) + s * model_->getNumCols(1),
							model_->getNumCols(1), solution_reco[augs_[s]]);

				/** collect solution */
				tssbdsub->solveRecourse(solution_, objval_reco, solution_reco, parNumCores_);
				for (int s = 0; s < model_->getNumScenarios(); ++s)
				{
					if (tssbdsub->excludedScenarios_[s]) continue;
					CoinCopyN(solution_reco[s], model_->getNumCols(1),
							solution_ + model_->getNumCols(0) + model_->getNumCols(1) * s);

				}

				/** compute primal bound */
				primalBound_ = 0.0;
				for (int j = 0; j < model_->getNumCols(0); ++j)
					primalBound_ += model_->getObjCore(0)[j] * solution_[j];
				for (int s = 0; s < model_->getNumScenarios(); ++s)
				{
					primalBound_ += objval_reco[s] * model_->getProbability()[s];
				}
				//printf("primalBound %e\n", primalBound_);

				/** free memory */
				FREE_ARRAY_PTR(objval_reco);
				FREE_2D_ARRAY_PTR(model_->getNumScenarios(), solution_reco);
			}
			break;
		}
		break;
	default:
		printf("Solution status (%d).\n", status_);
		break;
	}

	message_->print(1, " (%.2f sec)\n", CoinGetTimeOfDay() - tic);

	END_TRY_CATCH_RTN(FREE_MEMORY,STO_RTN_ERR)

	/** free memory */
	FREE_MEMORY

	return STO_RTN_OK;

#undef FREE_MEMORY
}

/** construct master problem */
STO_RTN_CODE TssBd::constructMasterProblem(TssBdSub * tssbdsub, double lowerbound)
{
#define FREE_MEMORY       \
	FREE_PTR(mat)         \
	FREE_ARRAY_PTR(clbd)  \
	FREE_ARRAY_PTR(cubd)  \
	FREE_ARRAY_PTR(ctype) \
	FREE_ARRAY_PTR(obj)   \
	FREE_ARRAY_PTR(rlbd)  \
	FREE_ARRAY_PTR(rubd)  \
	FREE_ARRAY_PTR(auxind)  \
	FREE_ARRAY_PTR(auxcoef)

	assert(model_);

	if (naux_ <= 0 || !obj_aux_ || !clbd_aux_ || !cubd_aux_)
	{
		printf("Warning: Auxiliary column information is required.\n");
		return 1;
	}

	/** master problem */
	CoinPackedMatrix * mat = NULL;
	double * clbd   = NULL;
	double * cubd   = NULL;
	double * obj    = NULL;
	char *   ctype  = NULL;
	double * rlbd   = NULL;
	double * rubd   = NULL;
	/** for recourse lower bound */
	int * auxind     = NULL;
	double * auxcoef = NULL;

	BGN_TRY_CATCH

	int ncols = model_->getNumCols(0) + naugs_ * model_->getNumCols(1) + naux_;

	/** number of integer variables in the core */
	int nIntegers = model_->getNumCoreIntegers();

	/** allocate memory */
	auxind = new int [ncols];
	auxcoef = new double [ncols];

	/** decompose model */
	STO_RTN_CHECK_THROW(model_->decompose(
			naugs_, augs_, naux_, clbd_aux_, cubd_aux_, obj_aux_,
			mat, clbd, cubd, ctype, obj, rlbd, rubd),
			"decompose", "TssModel");

	/** convert column types */
	if (parRelaxIntegrality_[0])
	{
		for (int j = 0; j < model_->getNumCols(0); ++j)
		{
			if (ctype[j] != 'C')
				nIntegers--;
			ctype[j] = 'C';
		}
	}
	if (parRelaxIntegrality_[1] && naugs_ > 0)
	{
		for (int j = 0; j < model_->getNumCols(1); ++j)
		{
			if (ctype[model_->getNumCols(0) + j] != 'C')
				nIntegers--;
		}
		CoinFillN(ctype + model_->getNumCols(0), naugs_ * model_->getNumCols(1), 'C');
	}

	if (nIntegers > 0)
	{
		si_ = new SolverInterfaceScip(par_);
		si_->setPrintLevel(CoinMin(parLogLevel_ + 1, 5));
	}
	else
	{
		si_ = new SolverInterfaceSpx(par_);
		si_->setPrintLevel(CoinMax(parLogLevel_ - 1, 0));
	}

	/** load problem data */
	si_->loadProblem(mat, clbd, cubd, obj, ctype, rlbd, rubd, "BdMaster");

	/** add row for lower bound */
	for (int j = 0; j < ncols; ++j)
		auxind[j] = j;
	CoinCopyN(obj, ncols, auxcoef);
	si_->addRow(ncols, auxind, auxcoef, lowerbound, COIN_DBL_MAX);

	/** save memory */
	FREE_MEMORY

	END_TRY_CATCH_RTN(FREE_MEMORY,STO_RTN_ERR)

	return STO_RTN_OK;

#undef FREE_MEMORY
}

/** configure Phase 1 */
STO_RTN_CODE TssBd::configureSLP()
{
	BGN_TRY_CATCH

	/** TODO: any configuration? */

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
}

/** to find a lower bound by solving a set of group subproblems */
STO_RTN_CODE TssBd::findLowerBound(double & lowerbound)
{
#define FREE_MEMORY       \
	FREE_PTR(mat)         \
	FREE_ARRAY_PTR(clbd)  \
	FREE_ARRAY_PTR(cubd)  \
	FREE_ARRAY_PTR(ctype) \
	FREE_ARRAY_PTR(obj)   \
	FREE_ARRAY_PTR(rlbd)  \
	FREE_ARRAY_PTR(rubd)

	/** problem */
	SolverInterface ** si = NULL;
	CoinPackedMatrix * mat = NULL;
	double * clbd = NULL;
	double * cubd = NULL;
	double * obj = NULL;
	char * ctype = NULL;
	double * rlbd = NULL;
	double * rubd = NULL;
	bool doNext = true;
	int augs[1];

	BGN_TRY_CATCH

	/** initialize lower bound */
	lowerbound = 0.0;

	/** create models */
	si = new SolverInterface * [model_->getNumScenarios()];
	for (int s = 0; s < model_->getNumScenarios(); ++s)
	{
		/** augmented scenario index */
		augs[0] = s;

		/** decompose model */
		STO_RTN_CHECK_THROW(
				model_->decompose(1, &augs[0], 0, NULL,
						NULL, NULL, mat, clbd, cubd, ctype, obj,
						rlbd, rubd), "decompose", "TssModel");
		//mat->verifyMtx(4);

		/** number of integer variables in the core */
		int nIntegers = model_->getNumCoreIntegers();

		/** relax second stage? */
		if (parRelaxIntegrality_[1])
		{
			/** relax integrality in the second stage */
			for (int j = 0; j < model_->getNumCols(1); ++j)
			{
				if (ctype[model_->getNumCols(0) + j] != 'C')
					nIntegers--;
				ctype[model_->getNumCols(0) + j] = 'C';
			}
		}
		//printf("ctype %s\n", ctype);

		/** adjust first-stage cost */
		for (int j = 0; j < model_->getNumCols(0); ++j)
			obj[j] *= model_->getProbability()[s]
;

		if (parInitLbAlgo_ == SEPARATE_MILP && nIntegers > 0)
		{
			si[s] = new SolverInterfaceScip(par_);
			//si[s]->setPrintLevel(CoinMax(CoinMin(parLogLevel_ - 1, 5),0));
			si[s]->setPrintLevel(0);
		}
		else
		{
			si[s] = new SolverInterfaceSpx(par_);
		}
		/** load problem */
		si[s]->loadProblem(mat, clbd, cubd, obj, ctype, rlbd, rubd, "BendersLowerBound");

		/** save memory */
		FREE_MEMORY
	}

	/** set number of cores to use */
#ifdef USE_OMP
	omp_set_num_threads(parNumCores_);
#pragma omp parallel for schedule(dynamic)
#endif
	for (int s = 0; s < model_->getNumScenarios(); ++s)
	{
		/** solve */
		si[s]->solve();
	}

	for (int s = 0; s < model_->getNumScenarios() && doNext == true; ++s)
	{
		/** solution status */
		status_ = si[s]->getStatus();

		/** get solution */
		switch (status_)
		{
		case STO_STAT_OPTIMAL:
		case STO_STAT_LIM_ITERorTIME:
			lowerbound += si[s]->getPrimalBound();
			break;
		case STO_STAT_DUAL_INFEASIBLE:
			/** need to use L-shaped method in this case */
			lowerbound = -COIN_DBL_MAX;
			doNext = false;
			break;
		default:
			doNext = false;
			break;
		}
	}

	END_TRY_CATCH_RTN(FREE_MEMORY;FREE_2D_PTR(model_->getNumScenarios(), si);,STO_RTN_ERR)

	FREE_2D_PTR(model_->getNumScenarios(), si);

	return STO_RTN_OK;
#undef FREE_MEMORY
}

/**
 * solve Phase 1
 *
 * TODO: Need more sophisticated cut management for efficiency.
 */
STO_RTN_CODE TssBd::solveSLP(TssBdSub * tssbdsub)
{
//#define CUT_AGGREGATION

	int iter = 0;
	bool repeat = true;
	OsiCuts cuts;
	int ncuts = 0;

	double objval = -COIN_DBL_MAX;
	double effectiveness   = 1E-8;
	double effectivenessLB = 1E-8;

	double stime_cpu  = 0.0;
	double stime_wall = 0.0;

	BGN_TRY_CATCH

	/** initial solve */
	si_->solve();

	int niters = 0;

	while (repeat)
	{
		/** solution status */
		if (si_->getStatus() == STO_STAT_OPTIMAL)
		{
			/** check improvement of objective value */
			effectiveness = fabs(objval - si_->getPrimalBound()) < 1E-10 ? effectiveness * 10 : effectivenessLB;

			/** optimal */
			objval = si_->getPrimalBound();

			/** print */
			if (parLogLevel_)
			{
				niters += si_->getIterationCount();
				message_->print(1,"Iteration %4d: objective function: %+E iteration %d\n",
						iter++, objval, niters);
			}
#ifdef TSSBENDERS_DEBUG
			PRINT_ARRAY(si_->getNumCols(), si_->getColSolution());
#endif

			/** stop at iteration limit */
			if (iter >= parIterLim_)
			{
				status_ = STO_STAT_STOPPED_ITER;
				repeat = false;
				break;
			}

			/** mark start times */
			stime_cpu  = CoinCpuTime();
			stime_wall = CoinGetTimeOfDay();

			/** generate cuts */
			tssbdsub->generateCuts(si_->getNumCols(), si_->getSolution(), &cuts);
			//cuts.printCuts();

			/** calculate effectiveness */
			for (int i = 0; i < cuts.sizeCuts(); ++i)
			{
				OsiRowCut * rc = cuts.rowCutPtr(i);
				if (!rc) continue;

				double violated = rc->violated(si_->getSolution());
				rc->setEffectiveness(violated);
			}

			/** mark elapsed times */
			stat_.cut_generation_time_cpu_phase1_  += CoinCpuTime() - stime_cpu;
			stat_.cut_generation_time_wall_phase1_ += CoinGetTimeOfDay() - stime_wall;

			/** add cuts */
			int nCutsAdded = si_->addCuts(cuts, effectiveness);

			/** move cuts */
			for (int i = 0; i < cuts.sizeCuts(); ++i)
			{
				OsiRowCut * rc = cuts.rowCutPtr(i);
				FREE_PTR(rc);
			}
			cuts.dumpCuts();

			/** if no cut added */
			ncuts += nCutsAdded;
			if (nCutsAdded == 0)
			{
				status_ = STO_STAT_OPTIMAL;
				repeat = false;
			}

			/** resolve master */
			si_->solve();
		}
		else
		{
			status_ = si_->getStatus();
			repeat = false;
		}
	}

	/** assign solution if no integer variable */
	if (status_ == STO_STAT_OPTIMAL ||
		status_ == STO_STAT_STOPPED_ITER)
	{
		/** primal bound */
		if (status_ == STO_STAT_OPTIMAL)
			primalBound_ = si_->getPrimalBound();
		else
			primalBound_ = COIN_DBL_MAX;

		/** dual bound */
		dualBound_ = si_->getPrimalBound();
	}

	/** statistics */
	numIterations_ = iter;

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
}

/** configure Phase 2 (e.g., Benders cut generator) */
STO_RTN_CODE TssBd::configureSMILP(TssBdSub * tssbdsub)
{
	BGN_TRY_CATCH

	/** retrieve solver interface for SCIP */
	SolverInterfaceScip * SiScip = dynamic_cast<SolverInterfaceScip*>(si_);
	assert(SiScip);

	/** get SCIP pointer */
	SCIP * scip = SiScip->getSCIP();

	/** create constraint handler */
	SCIPconshdlrBenders * conshdlr = new SCIPconshdlrBenders(scip, parCutPriority_);
	conshdlr->assignTssBdSub(tssbdsub);
	conshdlr->setOriginalVariables(SiScip->getNumCols(), SiScip->getSCIPvars());

	/** add constraint handler */
	SiScip->addConstraintHandler(conshdlr, true);

#if 0
	/** TODO: set branching priorities */
	if (model_->getPriorities() != NULL)
	{
	}

	/** TODO: set branching objects */
	if (model_->branchingHyperplanes_.size() > 0)
	{
	}
#endif

	/** set node limit */
	si_->setNodeLimit(parNodeLim_);

	/** set time limit */
	si_->setTimeLimit(time_remains_);

	/** set print level */
	si_->setPrintLevel(CoinMin(parLogLevel_ + 2, 5));

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
}

STO_RTN_CODE TssBd::solveSMILP()
{
	BGN_TRY_CATCH

	/** solve */
	si_->solve();

	/** solution status */
	status_ = si_->getStatus();

	switch (status_)
	{
	case STO_STAT_OPTIMAL:
	case STO_STAT_LIM_ITERorTIME:
	case STO_STAT_STOPPED_GAP:
	case STO_STAT_STOPPED_NODE:
	case STO_STAT_STOPPED_TIME:
		{
			primalBound_ = si_->getPrimalBound();
			dualBound_ = si_->getDualBound();
			break;
		}
		break;
	default:
		printf("Solution status (%d).\n", status_);
		break;
	}

	/** statistics */
	numIterations_ = si_->getIterationCount();
	numNodes_ = si_->getNumNodes();

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
}

/** set auxiliary variable data */
void TssBd::setAuxColData(int size, double * obj, double * clbd, double * cubd)
{
	FREE_ARRAY_PTR(obj_aux_)
	FREE_ARRAY_PTR(clbd_aux_)
	FREE_ARRAY_PTR(cubd_aux_)

	naux_ = size;
	obj_aux_  = new double [naux_];
	clbd_aux_ = new double [naux_];
	cubd_aux_ = new double [naux_];

	CoinCopyN(obj, naux_, obj_aux_);
	CoinCopyN(clbd, naux_, clbd_aux_);
	CoinCopyN(cubd, naux_, cubd_aux_);
}

/** set auxiliary variable data */
void TssBd::setAugScenarios(int size, int * indices)
{
	FREE_ARRAY_PTR(augs_)

	naugs_ = size;
	augs_  = new int [naugs_];

	CoinCopyN(indices, naugs_, augs_);
}

