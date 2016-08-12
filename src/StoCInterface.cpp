/*
 * StoCInterface.cpp
 *
 *  Created on: Oct 27, 2014
 *      Author: kibaekkim
 */

#include "StoCInterface.h"

#include "Utility/StoMacros.h"
#include "Solver/TssEval.h"
#include "Solver/TssDe.h"
#include "Solver/TssBd.h"
#include "Solver/TssDdMpi.h"

/* using __cplusplus */
#ifdef __cplusplus
extern "C" {
#endif

#define STO_API_CHECK_ENV(RTN)                \
	if (env == NULL) {                        \
		printf("Error: Null API pointer.\n"); \
		return RTN;                           \
	}
#define STO_API_CHECK_MODEL(RTN)                \
	STO_API_CHECK_ENV(RTN)                      \
	if (env->tss_ == NULL) {                    \
		printf("Error: Null model pointer.\n"); \
		return RTN;                             \
	}
#define STO_API_CHECK_SOLVER(RTN)                \
	STO_API_CHECK_ENV(RTN)                       \
	if (env->solver_ == NULL) {                  \
		printf("Error: Null solver pointer.\n"); \
		return RTN;                              \
	}

/** create API environment */
StoApiEnv * createEnv(void)
{
	return new StoApiEnv;
}

/** free API environment */
void freeEnv(StoApiEnv * env)
{
	FREE_PTR(env);
}

/** free new model */
void freeTssModel(StoApiEnv * env)
{
	STO_API_CHECK_ENV();
	FREE_PTR(env->tss_);
}

/** free solver */
void freeTssSolver(StoApiEnv * env)
{
	STO_API_CHECK_ENV();
	FREE_PTR(env->solver_);
}

/** get model pointer */
TssModel * getModelPtr(StoApiEnv * env)
{
	STO_API_CHECK_ENV(NULL);
	return env->tss_;
}

/** set number of scenarios */
void setNumberOfScenarios(
		StoApiEnv * env,  /**< pointer to API object */
		int         nscen /**< number of scenarios */)
{
	if (env->tss_ == NULL)
		env->tss_ = new TssModel;
	env->tss_->setNumberOfScenarios(nscen);
}

/** set dimensions */
void setDimensions(
		StoApiEnv * env,    /**< pointer to API object */
		const int   ncols1, /**< number of first-stage columns */
		const int   nrows1, /**< number of first-stage rows */
		const int   ncols2, /**< number of second-stage columns */
		const int   nrows2  /**< number of second-stage rows */)
{
	if (env->tss_ == NULL)
		env->tss_ = new TssModel;
	env->tss_->setDimensions(ncols1, nrows1, ncols2, nrows2);
}

/** read smps files */
void readSmps(StoApiEnv * env, const char * smps)
{
	if (env->tss_ == NULL)
		env->tss_ = new TssModel;
	env->tss_->readSmps(smps);
}

/** write the extensive form in Mps file */
void writeMps(StoApiEnv * env, const char * mps)
{
	STO_API_CHECK_MODEL();
	freeTssSolver(env);
	env->solver_ = new TssDe;
	env->solver_->loadModel(env->par_, env->tss_);
	dynamic_cast<TssDe*>(env->solver_)->writeMps(mps);
}

/** load first-stage problem */
void loadFirstStage(
		StoApiEnv *          env,   /**< pointer to API object */
		const CoinBigIndex * start, /**< start index for each row */
		const int *          index, /**< column indices */
		const double *       value, /**< constraint elements */
		const double *       clbd,  /**< column lower bounds */
		const double *       cubd,  /**< column upper bounds */
		const char *         ctype, /**< column types */
		const double *       obj,   /**< objective coefficients */
		const double *       rlbd,  /**< row lower bounds */
		const double *       rubd   /**< row upper bounds */)
{
	if (env->tss_ == NULL)
		env->tss_ = new TssModel;
	env->tss_->loadFirstStage(start, index, value, clbd, cubd, ctype, obj, rlbd, rubd);
}

/** load second-stage problem */
void loadSecondStage(
		StoApiEnv *          env,   /**< pointer to API object */
		const int            s,     /**< scenario index */
		const double         prob,  /**< probability */
		const CoinBigIndex * start, /**< start index for each row */
		const int *          index, /**< column indices */
		const double *       value, /**< constraint elements */
		const double *       clbd,  /**< column lower bounds */
		const double *       cubd,  /**< column upper bounds */
		const char *         ctype, /**< column types */
		const double *       obj,   /**< objective coefficients */
		const double *       rlbd,  /**< row lower bounds */
		const double *       rubd   /**< row upper bounds */)
{
	if (env->tss_ == NULL)
		env->tss_ = new TssModel;
	env->tss_->loadSecondStage(s, prob, start, index, value, clbd, cubd, ctype, obj, rlbd, rubd);
}


/** set branching priority */
void branchPriorities(
		StoApiEnv * env,       /**< pointer to API object */
		int *       priorities /**< A smaller number has more priority. */)
{
	if (env->tss_ == NULL)
		env->tss_ = new TssModel;
	env->tss_->setPriorities(priorities);
}

#if 0
/** add branching object */
void addBranchingObject(
		StoApiEnv * env,     /**< pointer to API object */
		int         nzcnt,   /**< number of nonzero elements */
		int *       indices, /**< indices */
		double *    values,  /**< hyper-plane coefficients */
		int         priority /**< branching priority */)
{
	if (env->tss_ == NULL)
		env->tss_ = new TssModel;
	env->tss_->addBranchingHyperplane(nzcnt, indices, values, priority);
}
#endif

/** evaluate solution */
void evaluateSolution(StoApiEnv * env, double * solution)
{
	STO_API_CHECK_MODEL();
	freeTssSolver(env);
	env->solver_ = new TssEval;
	env->solver_->loadModel(env->par_, env->tss_);
	TssEval * solver = dynamic_cast<TssEval*>(env->solver_);
	solver->setSolution(solution);
	env->solver_->solve();
}

/** solve deterministic equivalent model */
void solveDe(StoApiEnv * env)
{
	STO_API_CHECK_MODEL();
	freeTssSolver(env);
	env->solver_ = new TssDe;
	env->solver_->loadModel(env->par_, env->tss_);
	env->solver_->solve();
}

/** solve dual decomposition */
void solveDd(StoApiEnv * env, MPI_Comm comm)
{
	STO_API_CHECK_MODEL();
	freeTssSolver(env);
	env->solver_ = new TssDdMpi(comm, "");
	env->solver_->loadModel(env->par_, env->tss_);
	env->solver_->solve();
}

/** solve Benders decomposition */
void solveBd(
		StoApiEnv * env,     /**< pointer to API object */
		int         nauxvars /**< number of auxiliary variables (scenario clusters) */)
{
	STO_API_CHECK_MODEL();
	freeTssSolver(env);
	env->solver_ = new TssBd;
	env->solver_->loadModel(env->par_, env->tss_);

	/** set auxiliary variables */
	TssBd * tssbd = dynamic_cast<TssBd*>(env->solver_);
	if (tssbd == NULL)
	{
		printf("Error: TssSolver could not be converted to TssBd\n");
		return;
	}
	double * obj_aux = new double [nauxvars];
	double * clbd_aux = new double [nauxvars];
	double * cubd_aux = new double [nauxvars];
	CoinFillN(obj_aux, nauxvars, 1.0);
	CoinFillN(clbd_aux, nauxvars, -COIN_DBL_MAX);
	CoinFillN(cubd_aux, nauxvars, +COIN_DBL_MAX);
	tssbd->setAuxColData(nauxvars, obj_aux, clbd_aux, cubd_aux);

	/** set augmented scenarios */
	int numAugScenarios = env->par_->TssBdNumAugScenarios_;
	if (numAugScenarios > 0)
		tssbd->setAugScenarios(numAugScenarios, env->par_->TssBdAugScenarios_);

	/** relax second-stage integrality */
	setIntRelax(env, 1);

	FREE_ARRAY_PTR(obj_aux);
	FREE_ARRAY_PTR(clbd_aux);
	FREE_ARRAY_PTR(cubd_aux);

	env->solver_->solve();
}

/** set log level */
void setLogLevel(StoApiEnv * env, int level)
{
	STO_API_CHECK_ENV();
	env->par_->logLevel_ = level;
}

/** set number of cores */
void setNumCores(StoApiEnv * env, int num)
{
	STO_API_CHECK_ENV();
	env->par_->numCores_ = num;
}

/** set node limit */
void setNodeLimit(StoApiEnv * env, int num)
{
	STO_API_CHECK_ENV();
	env->par_->nodeLimit_ = num;
}

/** set iteration limit */
void setIterLimit(StoApiEnv * env, int num)
{
	STO_API_CHECK_ENV();
	env->par_->iterLimit_ = num;
}

/** set wallclock limit */
void setWallLimit(StoApiEnv * env, double lim)
{
	STO_API_CHECK_ENV();
	env->par_->wtimeLimit_ = lim;
}

/** set integrality relaxation */
void setIntRelax(StoApiEnv * env, int stage)
{
	STO_API_CHECK_ENV();
	env->par_->relaxIntegrality_[stage] = true;
}

/** set Benders augmented scenarios */
void setBdAugScenarios(StoApiEnv * env, int size, int * scenarios)
{
	STO_API_CHECK_ENV();
	env->par_->TssBdNumAugScenarios_ = size;
	env->par_->TssBdAugScenarios_ = new int [size];
	for (int i = 0; i < size; ++i)
		env->par_->TssBdAugScenarios_[i] = scenarios[i];
}

/** set Benders aggressiveness */
void setBendersAggressive(StoApiEnv * env, int aggressive)
{
	STO_API_CHECK_ENV();
	env->par_->TssBdBendersPriority_ = aggressive;
}

/** set a set of scenarios for the current process */
void setDdProcIdxSet(StoApiEnv * env, int size, int * scenarios)
{
	STO_API_CHECK_ENV();
	env->par_->TssDdNumProcIdx_ = size;
	env->par_->TssDdProcIdxSet_ = new int [size];
	for (int i = 0; i < size; ++i)
		env->par_->TssDdProcIdxSet_[i] = scenarios[i];
}

/** set parameter for adding feasibility cuts */
void setDdAddFeasCuts(StoApiEnv * env, int freq)
{
	STO_API_CHECK_ENV();
	env->par_->TssDdAddFeasCuts_ = freq;
}

/** set parameter for adding optimality cuts */
void setDdAddOptCuts(StoApiEnv * env, int freq)
{
	STO_API_CHECK_ENV();
	env->par_->TssDdAddOptCuts_ = freq;
}

/** set parameter for evaluating upper bound */
void setDdEvalUb(StoApiEnv * env, int freq)
{
	STO_API_CHECK_ENV();
	env->par_->TssDdEvalUb_ = freq;
}

/** set parameter for logging changes of distance of dual variables */
void setDdDualVarsLog(StoApiEnv * env, int yesNo)
{
	STO_API_CHECK_ENV();
	env->par_->TssDdDualVarsLog_ = yesNo;
}

/** set on/off DD recourse cache */
void setDdCacheRecourse(StoApiEnv * env, int yesNo)
{
	STO_API_CHECK_ENV();
	env->par_->TssDdCacheRecourse_ = yesNo;
}

/** set Lagrangian master solver */
void setDdMasterSolver(StoApiEnv * env, int type)
{
	STO_API_CHECK_ENV();
	env->par_->TssDdMasterSolver_ = type;
}

/** set DD stopping tolerance */
void setDdStoppingTolerance(StoApiEnv * env, double tol)
{
        STO_API_CHECK_ENV();
        env->par_->TssDdStoppingTol_ = tol;
}

/** set number of cuts per iteration added to master */
void setDdMasterNumCutsPerIter(StoApiEnv* env, int num)
{
	STO_API_CHECK_ENV();
	env->par_->TssDdMasterNumCutsPerIter_ = num;
}

/** set SCIP/display/freq */
void setScipDisplayFreq(StoApiEnv * env, int freq)
{
	STO_API_CHECK_ENV();
	env->par_->ScipDisplayFreq_ = freq;
}

/** set SCIP/limits/gap */
void setScipLimitsGap(StoApiEnv * env, double gap)
{
	STO_API_CHECK_ENV();
	env->par_->ScipLimitsGap_ = gap;
}

/** set SCIP/limits/time */
void setScipLimitsTime(StoApiEnv * env, double time)
{
	STO_API_CHECK_ENV();
	env->par_->ScipLimitsTime_ = time;
}

/**
 * Get functions
 */

/** get number of rows */
int getNumRows(StoApiEnv * env, int stage)
{
	STO_API_CHECK_MODEL(-1);
	return env->tss_->getNumRows(stage);
}

/** get number of columns */
int getNumCols(StoApiEnv * env, int stage)
{
	STO_API_CHECK_MODEL(-1);
	return env->tss_->getNumCols(stage);
}

/** get number of scenarios */
int getNumScenarios(StoApiEnv * env)
{
	STO_API_CHECK_MODEL(-1);
	return env->tss_->getNumScenarios();
}

/** get solution time */
double getSolutionTime(StoApiEnv * env)
{
	STO_API_CHECK_SOLVER(0.);
	return env->solver_->solutionTime_;
}

void getObjCoef(StoApiEnv * env, double * obj)
{
	CoinCopyN(env->tss_->getObjCore(0), env->tss_->getNumCols(0), obj);
	for (int s = 0; s < env->tss_->getNumScenarios(); ++s)
	{
		CoinCopyN(env->tss_->getObjCore(1), env->tss_->getNumCols(1), obj + env->tss_->getNumCols(0) + s * env->tss_->getNumCols(1));
		for (int j = 0; j < env->tss_->getObjScenario(s)->getNumElements(); ++j)
		{
			obj[s * env->tss_->getNumCols(1) + env->tss_->getObjScenario(s)->getIndices()[j]] = env->tss_->getObjScenario(s)->getElements()[j];
		}
		for (int j = 0; j < env->tss_->getNumCols(1); ++j)
			obj[env->tss_->getNumCols(0) + s * env->tss_->getNumCols(1) + j] *= env->tss_->getProbability()[s];
	}
}

/** get solution status */
int getSolutionStatus(StoApiEnv * env)
{
	STO_API_CHECK_SOLVER(0);
	return env->solver_->status_;
}

/** get objective value */
double getObjValue(StoApiEnv * env)
{
	STO_API_CHECK_SOLVER(0.0);
	return env->solver_->primalBound_;
}

/** get objective value */
double getPrimalBound(StoApiEnv * env)
{
	STO_API_CHECK_SOLVER(0.0);
	return env->solver_->primalBound_;
}

/** get objective value */
double getDualBound(StoApiEnv * env)
{
	STO_API_CHECK_SOLVER(0.0);
	return env->solver_->dualBound_;
}

/** get solution */
void getSolution(StoApiEnv * env, int num, double * solution)
{
	STO_API_CHECK_SOLVER();
	CoinCopyN(env->solver_->solution_, num, solution);
}

/** get number of iterations */
int getNumIterations(StoApiEnv * env)
{
	STO_API_CHECK_SOLVER(0);
	return env->solver_->numIterations_;
}

/** get number of nodes */
int getNumNodes(StoApiEnv * env)
{
	STO_API_CHECK_SOLVER(0);
	return env->solver_->numNodes_;
}

/** get number of infeasible solutions evaluated in DD */
int getDdNumInfeasSolutions(StoApiEnv * env)
{
	int num = 0;
	STO_API_CHECK_SOLVER(num);
	TssDdMpi * dd = dynamic_cast<TssDdMpi*>(env->solver_);
	if (dd)
		num = dd->nInfeasible_;
	return num;
}

/** get size of the iteration time array in DD */
int getDdIterTimeSize(StoApiEnv * env)
{
	int num = 0;
	STO_API_CHECK_SOLVER(num);
	TssDdMpi * dd = dynamic_cast<TssDdMpi*>(env->solver_);
	if (dd)
		num = dd->wtime_elapsed_.size();
	return num;
}

/** get solution time per iteration in DD */
void getDdIterTime(StoApiEnv * env, double * time)
{
	STO_API_CHECK_SOLVER();
	TssDdMpi * dd = dynamic_cast<TssDdMpi*>(env->solver_);
	if (dd)
	{
		for (unsigned int i = 0; i < dd->wtime_elapsed_.size(); ++i)
			time[i] = dd->wtime_elapsed_[i];
	}
}

/** get size of the master solution time array in DD */
int getDdMasterTimeSize(StoApiEnv * env)
{
	int num = 0;
	STO_API_CHECK_SOLVER(num);
	TssDdMpi * dd = dynamic_cast<TssDdMpi*>(env->solver_);
	if (dd)
		num = dd->wtime_master_.size();
	return num;
}

/** get master solution time per iteration in DD */
void getDdMasterTime(StoApiEnv * env, double * time)
{
	STO_API_CHECK_SOLVER();
	TssDdMpi * dd = dynamic_cast<TssDdMpi*>(env->solver_);
	if (dd)
	{
		for (unsigned int i = 0; i < dd->wtime_master_.size(); ++i)
			time[i] = dd->wtime_master_[i];
	}
}

/** get size of the subproblem solution time array in DD */
int getDdSubprobTimeSize(StoApiEnv * env)
{
	int num = 0;
	STO_API_CHECK_SOLVER(num);
	TssDdMpi * dd = dynamic_cast<TssDdMpi*>(env->solver_);
	if (dd)
		num = dd->wtime_subprob_.size();
	return num;
}

/** get subproblem solution time per iteration in DD */
void getDdSubprobTime(StoApiEnv * env, double * time)
{
	STO_API_CHECK_SOLVER();
	TssDdMpi * dd = dynamic_cast<TssDdMpi*>(env->solver_);
	if (dd)
	{
		for (unsigned int i = 0; i < dd->wtime_subprob_.size(); ++i)
			time[i] = dd->wtime_subprob_[i];
	}
}

/** get number of master objective values */
int getDdNumMasterObjValues(StoApiEnv * env)
{
	STO_API_CHECK_SOLVER(0);
	TssDdMpi * dd = dynamic_cast<TssDdMpi*>(env->solver_);
	if (dd)
	{
		return dd->objval_master_.size();
	}
	return 0;
}

/** get number of subproblem objective values */
int getDdNumSubproblemObjValues(StoApiEnv * env)
{
	STO_API_CHECK_SOLVER(0);
	TssDdMpi * dd = dynamic_cast<TssDdMpi*>(env->solver_);
	if (dd)
	{
		return dd->objval_subprob_.size();
	}
	return 0;
}

/** get number of primal bounds */
int getDdNumPrimalBounds(StoApiEnv * env)
{
	STO_API_CHECK_SOLVER(0);
	TssDdMpi * dd = dynamic_cast<TssDdMpi*>(env->solver_);
	if (dd)
	{
		return dd->primalBounds_.size();
	}
	return 0;
}

/** get number of dual bounds */
int getDdNumDualBounds(StoApiEnv * env)
{
	STO_API_CHECK_SOLVER(0);
	TssDdMpi * dd = dynamic_cast<TssDdMpi*>(env->solver_);
	if (dd)
	{
		return dd->dualBounds_.size();
	}
	return 0;
}

/** get history of master objective values */
void getDdMasterObjValues(StoApiEnv * env, double * vals)
{
	STO_API_CHECK_SOLVER();
	TssDdMpi * dd = dynamic_cast<TssDdMpi*>(env->solver_);
	if (dd)
	{
		for (unsigned int i = 0; i < dd->objval_master_.size(); ++i)
			vals[i] = dd->objval_master_[i];
	}
}

/** get history of subproblem objective values */
void getDdSubproblemObjValues(StoApiEnv * env, double * vals)
{
	STO_API_CHECK_SOLVER();
	TssDdMpi * dd = dynamic_cast<TssDdMpi*>(env->solver_);
	if (dd)
	{
		for (unsigned int i = 0; i < dd->objval_subprob_.size(); ++i)
			vals[i] = dd->objval_subprob_[i];
	}
}

/** get history of primal bounds */
void getDdPrimalBounds(StoApiEnv * env, double * vals)
{
	STO_API_CHECK_SOLVER();
	TssDdMpi * dd = dynamic_cast<TssDdMpi*>(env->solver_);
	if (dd)
	{
		for (unsigned int i = 0; i < dd->primalBounds_.size(); ++i)
			vals[i] = dd->primalBounds_[i];
	}
}

/** get history of dual bounds */
void getDdDualBounds(StoApiEnv * env, double * vals)
{
	STO_API_CHECK_SOLVER();
	TssDdMpi * dd = dynamic_cast<TssDdMpi*>(env->solver_);
	if (dd)
	{
		for (unsigned int i = 0; i < dd->dualBounds_.size(); ++i)
			vals[i] = dd->dualBounds_[i];
	}
}

/** get total cpu time */
double getDdCpuTime(StoApiEnv * env)
{
	STO_API_CHECK_SOLVER(0.);
	TssDdMpi * dd = dynamic_cast<TssDdMpi*>(env->solver_);
	if (dd)
	{
		return dd->ctime_start_;
	}
	return 0;
}

int getDdNumChangesOfMultiplier(StoApiEnv * env)
{
	STO_API_CHECK_SOLVER(0);
	TssDdMpi * dd = dynamic_cast<TssDdMpi*>(env->solver_);
	if (dd)
	{
		return dd->changesOfMultiplier_.size();
	}
	return 0;
}

/** get changes of multiplier */
void getDdChangesOfMultiplier(StoApiEnv * env, double * changes)
{
	STO_API_CHECK_SOLVER();
	TssDdMpi * dd = dynamic_cast<TssDdMpi*>(env->solver_);
	if (dd)
	{
		for (unsigned int i = 0; i < dd->changesOfMultiplier_.size(); ++i)
			changes[i] = dd->changesOfMultiplier_[i];
	}
}

/** get DD master time */
double getDdMasterTotalTime(StoApiEnv * env)
{
	STO_API_CHECK_SOLVER(0.);
	TssDdMpi * dd = dynamic_cast<TssDdMpi*>(env->solver_);
	if (dd)
	{
		return dd->total_time_master_;
	}
	return 0;
}

/** get DD lower bounding time */
double getDdLbTotalTime(StoApiEnv * env)
{
	STO_API_CHECK_SOLVER(0.);
	TssDdMpi * dd = dynamic_cast<TssDdMpi*>(env->solver_);
	if (dd)
	{
		return dd->total_time_lb_;
	}
	return 0;
}

/** get DD upper bounding time */
double getDdUbTotalTime(StoApiEnv * env)
{
	STO_API_CHECK_SOLVER(0.);
	TssDdMpi * dd = dynamic_cast<TssDdMpi*>(env->solver_);
	if (dd)
	{
		return dd->total_time_ub_;
	}
	return 0;
}

/** get DD cut generation time */
double getDdCgTotalTime(StoApiEnv * env)
{
	STO_API_CHECK_SOLVER(0.);
	TssDdMpi * dd = dynamic_cast<TssDdMpi*>(env->solver_);
	if (dd)
	{
		return dd->total_time_cg_;
	}
	return 0;
}

/**
 * Misc
 */

/** print model */
void printModel(StoApiEnv * env)
{
	STO_API_CHECK_MODEL();
	env->tss_->__printData();
}

#undef STO_API_CHECK_MODEL
#undef STO_API_CHECK_ENV

#ifdef __cplusplus
}
#endif


