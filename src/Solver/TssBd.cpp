/*
 * TssBd.cpp
 *
 *  Created on: Dec 5, 2014
 *      Author: kibaekkim
 */

//#define DSP_DEBUG

/** DSP */
#include "Solver/TssBd.h"
#include "Solver/SCIPconshdlrBenders.h"
#include "Solver/SolverInterfaceSpx.h"
#include "Solver/SolverInterfaceClp.h"
#include "Solver/SolverInterfaceScip.h"
#include "Utility/StoMessage.h"

TssBd::TssBd() :
	TssSolver(),
	si_(NULL),
	naugs_(0),
	augs_(NULL),
	naux_(0),
	obj_aux_(NULL),
	clbd_aux_(NULL),
	cubd_aux_(NULL)
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

/** solve */
STO_RTN_CODE TssBd::solve()
{
#define FREE_MEMORY    \
	FREE_PTR(tssbdsub)

	assert(model_);

	/** subproblems */
	const double * probability = NULL; /**< probability */
	TssBdSub *     tssbdsub    = NULL; /**< cut generator */
	double lowerbound = 0.0;

	BGN_TRY_CATCH

	/** solution time */
	double swtime = CoinGetTimeOfDay();

	/** time limit */
	time_remains_ = par_->wtimeLimit_;

	/** get number of scenarios */
	probability = model_->getProbability();

	/** configure Benders cut generator */
	/** This does NOT make deep copies. So, DO NOT RELEASE POINTERS. */
	tssbdsub = new TssBdSub(par_);
	for (int s = 0; s < model_->getNumScenarios(); ++s)
		tssbdsub->scenarios_.push_back(s);
	STO_RTN_CHECK_THROW(tssbdsub->loadProblem(model_, naugs_, augs_, naux_), "loadProblem", "TssBdSub");

	if (par_->logLevel_ > 0) printf("Phase 1:\n");

	/** solution time */
	double stime = clockType_ ? CoinGetTimeOfDay() : CoinCpuTime();

	/** find lower bound */
	/** TODO: replace this with TssDd */
	STO_RTN_CHECK_THROW(findLowerBound(probability, lowerbound), "findLowerBound", "TssBd");

	/** We should have a lower bound here. */
	if (par_->logLevel_ > 0) printf(" -> Found lower bound %e\n", lowerbound);

	/** construct master problem */
	STO_RTN_CHECK_THROW(constructMasterProblem(tssbdsub, lowerbound), "constructMasterProblem", "TssBd");

	if (par_->logLevel_ > 0) printf("Phase 2:\n");

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
				tssbdsub->solveRecourse(solution_, objval_reco, solution_reco, par_->numCores_);
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
	if (par_->relaxIntegrality_[0])
	{
		for (int j = 0; j < model_->getNumCols(0); ++j)
		{
			if (ctype[j] != 'C')
				nIntegers--;
			ctype[j] = 'C';
		}
	}
	if (par_->relaxIntegrality_[1] && naugs_ > 0)
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
		si_->setPrintLevel(CoinMin(par_->logLevel_ + 1, 5));
	}
	else
	{
		si_ = new SolverInterfaceSpx(par_);
		si_->setPrintLevel(CoinMax(par_->logLevel_ - 1, 0));
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
STO_RTN_CODE TssBd::findLowerBound(
		const double * probability,
		double & lowerbound)
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
		if (par_->relaxIntegrality_[1])
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
			obj[j] *= probability[s];
#if 1
		if (nIntegers > 0)
		{
			si[s] = new SolverInterfaceScip(par_);
			si[s]->setPrintLevel(par_->logLevel_ + 1);
		}
		else
#endif
		{
			si[s] = new SolverInterfaceSpx(par_);
			si[s]->setPrintLevel(CoinMax(par_->logLevel_ - 1, 0));
		}
		/** load problem */
		si[s]->loadProblem(mat, clbd, cubd, obj, ctype, rlbd, rubd, "BendersLowerBound");

		/** save memory */
		FREE_MEMORY
	}

	/** set number of cores to use */
#ifdef USE_OMP
	omp_set_num_threads(par_->numCores_);
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
			if (par_->logLevel_)
			{
				niters += si_->getIterationCount();
				printf("Iteration %4d: objective function: %+E iteration %d\n",
						iter++, objval, niters);
			}
#ifdef TSSBENDERS_DEBUG
			PRINT_ARRAY(si_->getNumCols(), si_->getColSolution());
#endif

			/** stop at iteration limit */
			if (iter >= par_->iterLimit_)
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
	SCIPconshdlrBenders * conshdlr = new SCIPconshdlrBenders(scip, par_->TssBdBendersPriority_);
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
	si_->setNodeLimit(par_->nodeLimit_);

	/** set time limit */
	si_->setTimeLimit(time_remains_);

	/** set print level */
	si_->setPrintLevel(CoinMin(par_->logLevel_ + 2, 5));

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

#if 0
/** construct subproblems */
STO_RTN_CODE TssBd::constructSubProblem(
		double * probability,
		CoinPackedMatrix **&   mat_tech,
		OsiSolverInterface **& si_reco)
{
#define FREE_MEMORY                 \
	FREE_PTR(mat_tech_l)            \
	FREE_PTR(mat_reco)              \
	FREE_ARRAY_PTR(clbd_reco)       \
	FREE_ARRAY_PTR(cubd_reco)       \
	FREE_ARRAY_PTR(ctype_reco)      \
	FREE_ARRAY_PTR(obj_reco)        \
	FREE_ARRAY_PTR(rlbd_reco)       \
	FREE_ARRAY_PTR(rubd_reco)       \
	FREE_ARRAY_PTR(intind_reco)     \
	FREE_ARRAY_PTR(shifted_indices)

	assert(model_);

	int i, s;
	std::vector<int> incscn, excscn; /** scenario indices to be included and excluded */

	/** subproblems */
	CoinPackedMatrix * mat_tech_l = NULL;
	CoinPackedMatrix * mat_reco = NULL;
	double * clbd_reco   = NULL;
	double * cubd_reco   = NULL;
	double * obj_reco    = NULL;
	char *   ctype_reco  = NULL;
	double * rlbd_reco   = NULL;
	double * rubd_reco   = NULL;
	int      intlen_reco = 0;
	int *    intind_reco = NULL;
	int *    shifted_indices = NULL;

	BGN_TRY_CATCH

	/** get number of scenarios */
	int nscen = model_->getNumScenarios();

	/** get incscn and excscn */
	for (s = 0, i = 0; s < nscen; ++s)
	{
		/** for scenario augmentation */
		if (naugs_ > 0 && i < naugs_)
		{
			while (i < naugs_ && augs_[i] < s)
			{
				i++;
			}
			if (i < naugs_ && augs_[i] == s)
			{
				/** put in the set of excluded scenarios */
				excscn.push_back(s);
				continue;
			}
		}

		/** put in the set of included scenarios */
		incscn.push_back(s);
	}

	/** for included scenarios */
	for (unsigned int s = 0; s < incscn.size(); ++s)
	{
		/** get probability */
		probability[s] = model_->getProbability()[incscn[s]];

		/** copy recourse problem */
		STO_RTN_CHECK_THROW(model_->copyRecoProb(incscn[s], mat_tech[s],
				mat_reco, clbd_reco, cubd_reco, ctype_reco,
				obj_reco, rlbd_reco, rubd_reco), "copyRecoProb", "TssModel")

		/** column types */
		if (relax_[1]) CoinFillN(ctype_reco, model_->getNumCols(1), 'C');
		convertColTypes(mat_reco->getNumCols(), ctype_reco, intlen_reco, intind_reco);

		/** creating solver interface */
		si_reco[s] = new OsiClpSolverInterface;
		si_reco[s]->loadProblem(*mat_reco, clbd_reco, cubd_reco, obj_reco, rlbd_reco, rubd_reco);
		si_reco[s]->setInteger(intind_reco, intlen_reco);
		si_reco[s]->messageHandler()->setLogLevel(0);

		//printf("s %d mat_reco %p si_reco %p\n", s, mat_tech[s], si_reco[s]);

		/** free memory */
		FREE_MEMORY
	}

#if 0
	/** for excluded scenarios, create augmented rows */
	for (unsigned int s = 0; s < excscn.size(); ++s)
	{
		printf("Create augmented row for scenario %d\n", excscn[s]);
		/** copy recourse problem */
		status = model_->copyRecoProb(excscn[s], mat_tech_l, mat_reco, clbd_reco, cubd_reco, ctype_reco,
				obj_reco, rlbd_reco, rubd_reco);
		if (status) CoinError("Failed to copy recourse problem data", "constructSubProblem", "TssBenders");

		/** merge */
		mat_tech_l->rightAppendPackedMatrix(*mat_reco);
		mat_tech_l->verifyMtx(4);

		/** allocate temporary memory */
		shifted_indices = new int [mat_tech_l->getNumElements()];

		/** assign indices */
		CoinCopyN(mat_tech_l->getIndices(), mat_tech_l->getNumElements(), shifted_indices);

		/** shift indices */
		if (s > 0)
		{
			for (int i = 0; i < mat_tech_l->getNumRows(); ++i)
			{
				for (int j = mat_tech_l->getVectorStarts()[i]; j < mat_tech_l->getVectorStarts()[i+1]; ++j)
				{
					model_->shiftVecIndices(mat_tech_l->getVectorSize(i),
							shifted_indices + mat_tech_l->getVectorFirst(i),
							model_->getNumCols(1) * s,
							model_->getNumCols(0));
				}
			}
		}

		/** add rows */
		si_->addRows(mat_tech_l->getNumRows(), mat_tech_l->getVectorStarts(),
				shifted_indices, mat_tech_l->getElements(), rlbd_reco, rubd_reco);

		/** free memory */
		FREE_MEMORY
	}
#endif

//	/** set subproblem info to LsSolverInterface */
//	si_->loadSubproblems(
//			nscen - naugs_, probability,
//			const_cast<const CoinPackedMatrix **>(mat_tech),
//			const_cast<const OsiSolverInterface **>(si_reco),
//			naux_);

	END_TRY_CATCH_RTN(FREE_MEMORY,STO_RTN_ERR)

	return STO_RTN_OK;

#undef FREE_MEMORY
}
#endif

