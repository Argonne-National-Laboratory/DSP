/*
 * StoCInterface.cpp
 *
 *  Created on: Oct 27, 2014
 *      Author: kibaekkim
 */

#include "StoCInterface.h"

#include "Utility/StoMacros.h"
#include "Solver/TssEval.h"
#include "Solver/TssBd.h"
#include "Solver/DecDdMpi.h"
#include "Solver/DecDe.h"
#include "Solver/DecTssSolver.h"
#include "Model/DecTssModel.h"
#include "Model/DecDetModel.h"


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
	if (env->model_ == NULL) {                    \
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
void freeModel(StoApiEnv * env)
{
	STO_API_CHECK_ENV();
	FREE_PTR(env->model_);
}

/** free solver */
void freeSolver(StoApiEnv * env)
{
	STO_API_CHECK_ENV();
	FREE_PTR(env->solver_);
}

/** If current model is stochastic, return the model as a TssModel object. If no model exists, create one. */
TssModel * getTssModel(StoApiEnv * env)
{
	if (env->model_ == NULL)
		env->model_ = new DecTssModel;
	/** TODO: Should fail gracefully */
	if (env->model_->isStochastic())
	{
		TssModel * tss;
		try
		{
			tss = dynamic_cast<TssModel *>(env->model_);
		}
		catch (const std::bad_cast& e)
		{
			printf("Error: Model claims to be stochastic when it is not");
			return NULL;
		}
		return tss;
	}
	else
	{
		printf("Error: Attempted to access feature only supported by stochastic models with a general decomposition model\n");
		return NULL;
	}
}

/** prepare decomposition model to be solved; returns false if there is an error */
bool prepareDecModel(StoApiEnv * env)
{
	if (!env->model_->isStochastic() && env->decdata_ == NULL)
	{
		printf("Error: General decomposition models must be accompanied by coupling constraints\n");
		return false;
	}

	if (env->decdata_ != NULL)
	{
		if (env->model_->isStochastic())
		{
			printf("Decomposition data for a stochastic model supplied: converting model to extensive form\n");
			DetModel * det;
			STO_RTN_CHECK_THROW(getTssModel(env)->copyDeterministicEquivalent(det),
					"copyDeterministicEquivalent", "TssModel");
			env->model_ = new DecDetModel(det, env->decdata_);
			/** let default settings for subproblems take place */
			env->par_->TssDdNumProcIdx_ = 0;
			env->par_->TssDdProcIdxSet_ = NULL;
		}
		else
		{
			/** deterministic case with decdata: update model with decdata */
			DecDetModel * decDet;
			try
			{
				decDet = dynamic_cast<DecDetModel *>(env->model_);
			}
			catch (const std::bad_cast& e)
			{
				printf("Error: Deterministic model not of proper type");
				return false;
			}
			decDet->setDecData(env->decdata_);
		}

		if (env->par_->relaxIntegrality_[0] || env->par_->relaxIntegrality_[1])
		{
			printf("Warning: Relaxing stage integrality only supported in stochastic model; feature disabled\n");
		}
	}
	return true;
}

/** get model pointer */
DecModel * getModelPtr(StoApiEnv * env)
{
	STO_API_CHECK_ENV(NULL);
	return env->model_;
}

/** set number of scenarios */
void setNumberOfScenarios(
		StoApiEnv * env,  /**< pointer to API object */
		int         nscen /**< number of scenarios */)
{
	getTssModel(env)->setNumberOfScenarios(nscen);
}

/** set dimensions */
void setDimensions(
		StoApiEnv * env,    /**< pointer to API object */
		const int   ncols1, /**< number of first-stage columns */
		const int   nrows1, /**< number of first-stage rows */
		const int   ncols2, /**< number of second-stage columns */
		const int   nrows2  /**< number of second-stage rows */)
{
	getTssModel(env)->setDimensions(ncols1, nrows1, ncols2, nrows2);
}

/** read smps files */
void readSmps(StoApiEnv * env, const char * smps)
{
	getTssModel(env)->readSmps(smps);
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
	getTssModel(env)->loadFirstStage(start, index, value, clbd, cubd, ctype, obj, rlbd, rubd);
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
	getTssModel(env)->loadSecondStage(s, prob, start, index, value, clbd, cubd, ctype, obj, rlbd, rubd);
}

/** load deterministic problem */
void loadDeterministic(
		StoApiEnv *          env,    /**< pointer to API object */
		const CoinBigIndex * start,  /**< start index for each row */
		const int *          index,  /**< column indices */
		const double *       value,  /**< constraint elements */
		const int            numels, /**< number of elements in index and value */
		const int            ncols,  /**< number of columns */
		const int            nrows,  /**< number of rows */
		const double *       clbd,   /**< column lower bounds */
		const double *       cubd,   /**< column upper bounds */
		const char *         ctype,  /**< column types */
		const double *       obj,    /**< objective coefficients */
		const double *       rlbd,   /**< row lower bounds */
		const double *       rubd    /**< row upper bounds */)
{
	if (env->model_ != NULL)
		printf("Warning: Replacing an already loaded model\n");
	DetModel * det = new DetModel(start, index, value, numels, ncols, nrows, clbd, cubd, ctype, obj, rlbd, rubd);
	env->model_ = new DecDetModel(det, NULL);
}

/** load parameters for a custom decomposition of the problem */
void loadDecomposition(
		StoApiEnv * env,            /**< pointer to API object */
		int         nsubprobs,      /**< number of subproblems */
		int         ncols,          /**< number of columns */
		int         ncoupling,      /**< number of coupling constraints */
		int *       varPartition,   /**< partition of columns into subproblems */
		int *       couplingStarts, /**< indices in cols at which each coupling constraint starts */
		int *       couplingCols,   /**< variables of each coupling constraint left-hand side */
		double *    couplingCoeffs, /**< coefficients of each coupling constraint left-hand side */
		char *      couplingSenses, /**< senses of each coupling constraint */
		double *    couplingRhs     /**< right-hand sides of each coupling constraint */)
{
	if (env->model_ == NULL)
	{
		printf("Error: Model needs to be loaded before specifying decomposition parameters.\n");
		return;
	}

	env->decdata_ = new DecData(nsubprobs, ncols, ncoupling, varPartition, couplingStarts,
		couplingCols, couplingCoeffs, couplingSenses, couplingRhs);
}


/** set branching priority */
void branchPriorities(
		StoApiEnv * env,       /**< pointer to API object */
		int *       priorities /**< A smaller number has more priority. */)
{
	getTssModel(env)->setPriorities(priorities);
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
	getTssModel(env)->addBranchingHyperplane(nzcnt, indices, values, priority);
}
#endif

/** evaluate solution */
void evaluateSolution(StoApiEnv * env, double * solution)
{
	STO_API_CHECK_MODEL();
	freeSolver(env);

	if (!env->model_->isStochastic())
	{
		printf("Error: Solution evaluation is currently only supported by stochastic models\n");
		return;
	}

	env->solver_ = new DecTssSolver(new TssEval);
	env->solver_->loadModel(env->par_, new DecTssModel(*getTssModel(env)));
	TssEval * solver = dynamic_cast<TssEval*>(env->solver_);
	solver->setSolution(solution);
	env->solver_->solve();
}

/** solve deterministic equivalent model */
void solveDe(StoApiEnv * env)
{
	STO_API_CHECK_MODEL();
	freeSolver(env);

	if (!prepareDecModel(env))
		return;

	env->solver_ = new DecDe;
	env->solver_->loadModel(env->par_, env->model_);
	env->solver_->solve();
}

/** solve dual decomposition */
void solveDd(StoApiEnv * env, MPI_Comm comm)
{
	STO_API_CHECK_MODEL();
	freeSolver(env);

	if (!prepareDecModel(env))
		return;

	env->solver_ = new DecDdMpi(comm, "");
	env->solver_->loadModel(env->par_, env->model_);
	env->solver_->solve();
}

/** solve Benders decomposition */
void solveBd(
		StoApiEnv * env,     /**< pointer to API object */
		int         nauxvars /**< number of auxiliary variables (scenario clusters) */)
{
	STO_API_CHECK_MODEL();
	freeSolver(env);

	if (!env->model_->isStochastic())
	{
		printf("Error: Benders decomposition is currently only supported by stochastic models\n");
		return;
	}

	TssBd * tssbd = new TssBd;
	env->solver_ = new DecTssSolver(tssbd);
	env->solver_->loadModel(env->par_, new DecTssModel(*getTssModel(env)));

	/** set auxiliary variables */
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

/** set trust region size */
void setDdTrustRegionSize(StoApiEnv* env, int num)
{
	STO_API_CHECK_ENV();
	env->par_->TssDdTrustRegionSize_ = num;
}

/** set whether trust region decrease should be disabled */
void setDdDisableTrustRegionDecrease(StoApiEnv* env, bool yesNo)
{
	STO_API_CHECK_ENV();
	env->par_->TssDdDisableTrustRegionDecrease_ = yesNo;
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
	return getTssModel(env)->getNumRows(stage);
}

/** get number of columns */
int getNumCols(StoApiEnv * env, int stage)
{
	STO_API_CHECK_MODEL(-1);
	return getTssModel(env)->getNumCols(stage);
}

/** get number of columns */
int getTotalNumCols(StoApiEnv * env)
{
	STO_API_CHECK_MODEL(-1);
	return env->model_->getFullModelNumCols();
}

/** get number of scenarios */
int getNumScenarios(StoApiEnv * env)
{
	STO_API_CHECK_MODEL(-1);
	return getTssModel(env)->getNumScenarios();
}

/** get number of subproblems */
int getNumSubproblems(StoApiEnv * env)
{
	STO_API_CHECK_MODEL(-1);
	return env->model_->getNumSubproblems();
}

/** get solution time */
double getSolutionTime(StoApiEnv * env)
{
	STO_API_CHECK_SOLVER(0.);
	return env->solver_->solutionTime_;
}

void getObjCoef(StoApiEnv * env, double * obj)
{
	STO_API_CHECK_SOLVER();
	return env->model_->getObjCoef(obj);
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
	DecDdMpi * dd = dynamic_cast<DecDdMpi*>(env->solver_);
	if (dd)
		num = dd->nInfeasible_;
	return num;
}

/** get size of the iteration time array in DD */
int getDdIterTimeSize(StoApiEnv * env)
{
	int num = 0;
	STO_API_CHECK_SOLVER(num);
	DecDdMpi * dd = dynamic_cast<DecDdMpi*>(env->solver_);
	if (dd)
		num = dd->wtime_elapsed_.size();
	return num;
}

/** get solution time per iteration in DD */
void getDdIterTime(StoApiEnv * env, double * time)
{
	STO_API_CHECK_SOLVER();
	DecDdMpi * dd = dynamic_cast<DecDdMpi*>(env->solver_);
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
	DecDdMpi * dd = dynamic_cast<DecDdMpi*>(env->solver_);
	if (dd)
		num = dd->wtime_master_.size();
	return num;
}

/** get master solution time per iteration in DD */
void getDdMasterTime(StoApiEnv * env, double * time)
{
	STO_API_CHECK_SOLVER();
	DecDdMpi * dd = dynamic_cast<DecDdMpi*>(env->solver_);
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
	DecDdMpi * dd = dynamic_cast<DecDdMpi*>(env->solver_);
	if (dd)
		num = dd->wtime_subprob_.size();
	return num;
}

/** get subproblem solution time per iteration in DD */
void getDdSubprobTime(StoApiEnv * env, double * time)
{
	STO_API_CHECK_SOLVER();
	DecDdMpi * dd = dynamic_cast<DecDdMpi*>(env->solver_);
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
	DecDdMpi * dd = dynamic_cast<DecDdMpi*>(env->solver_);
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
	DecDdMpi * dd = dynamic_cast<DecDdMpi*>(env->solver_);
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
	DecDdMpi * dd = dynamic_cast<DecDdMpi*>(env->solver_);
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
	DecDdMpi * dd = dynamic_cast<DecDdMpi*>(env->solver_);
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
	DecDdMpi * dd = dynamic_cast<DecDdMpi*>(env->solver_);
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
	DecDdMpi * dd = dynamic_cast<DecDdMpi*>(env->solver_);
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
	DecDdMpi * dd = dynamic_cast<DecDdMpi*>(env->solver_);
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
	DecDdMpi * dd = dynamic_cast<DecDdMpi*>(env->solver_);
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
	DecDdMpi * dd = dynamic_cast<DecDdMpi*>(env->solver_);
	if (dd)
	{
		return dd->ctime_start_;
	}
	return 0;
}

int getDdNumChangesOfMultiplier(StoApiEnv * env)
{
	STO_API_CHECK_SOLVER(0);
	DecDdMpi * dd = dynamic_cast<DecDdMpi*>(env->solver_);
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
	DecDdMpi * dd = dynamic_cast<DecDdMpi*>(env->solver_);
	if (dd)
	{
		for (unsigned int i = 0; i < dd->changesOfMultiplier_.size(); ++i)
			changes[i] = dd->changesOfMultiplier_[i];
	}
}

/**
 * Misc
 */

/** print model */
void printModel(StoApiEnv * env)
{
	STO_API_CHECK_MODEL();
	env->model_->__printData();
}

#undef STO_API_CHECK_MODEL
#undef STO_API_CHECK_ENV

#ifdef __cplusplus
}
#endif


