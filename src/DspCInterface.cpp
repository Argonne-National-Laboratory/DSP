/*
 * DspCInterface.cpp
 *
 *  Created on: Oct 27, 2014
 *      Author: kibaekkim
 */

#include "DspCInterface.h"
#include "Utility/StoMacros.h"
#include "Solver/Deterministic/DeDriver.h"
#include "Solver/Benders/BdDriver.h"
#include "Solver/DualDecomp/DdDriver.h"
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
DspApiEnv * createEnv(void)
{
	return new DspApiEnv;
}

/** free API environment */
void freeEnv(DspApiEnv * env)
{
	FREE_PTR(env);
}

/** free new model */
void freeModel(DspApiEnv * env)
{
	STO_API_CHECK_ENV();
	FREE_PTR(env->model_);
}

/** free solver */
void freeSolver(DspApiEnv * env)
{
	STO_API_CHECK_ENV();
	FREE_PTR(env->solver_);
}

/** If current model is stochastic, return the model as a TssModel object. If no model exists, create one. */
TssModel * getTssModel(DspApiEnv * env)
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
bool prepareDecModel(DspApiEnv * env)
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

		if (env->par_->getBoolPtrParam("RELAX_INTEGRALITY")[0] ||
				env->par_->getBoolPtrParam("RELAX_INTEGRALITY")[1])
		{
			printf("Warning: Relaxing stage integrality only supported in stochastic model; feature disabled\n");
		}
	}
	return true;
}

/** get model pointer */
DecModel * getModelPtr(DspApiEnv * env)
{
	STO_API_CHECK_ENV(NULL);
	return env->model_;
}

/** set number of scenarios */
void setNumberOfScenarios(
		DspApiEnv * env,  /**< pointer to API object */
		int         nscen /**< number of scenarios */)
{
	getTssModel(env)->setNumberOfScenarios(nscen);
}

/** set dimensions */
void setDimensions(
		DspApiEnv * env,    /**< pointer to API object */
		const int   ncols1, /**< number of first-stage columns */
		const int   nrows1, /**< number of first-stage rows */
		const int   ncols2, /**< number of second-stage columns */
		const int   nrows2  /**< number of second-stage rows */)
{
	getTssModel(env)->setDimensions(ncols1, nrows1, ncols2, nrows2);
}

/** read smps files */
void readSmps(DspApiEnv * env, const char * smps)
{
	getTssModel(env)->readSmps(smps);
}

/** load first-stage problem */
void loadFirstStage(
		DspApiEnv *          env,   /**< pointer to API object */
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
		DspApiEnv *          env,   /**< pointer to API object */
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
		DspApiEnv *          env,    /**< pointer to API object */
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
		DspApiEnv * env,            /**< pointer to API object */
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

/** set initial solutions
 * This function can be called multiple times for multiple initial solutions. */
void setSolution(
		DspApiEnv * env,     /**< pointer to API object */
		int         size,    /**< size of solution array */
		double *    solution /**< solution to set */)
{
	getTssModel(env)->setSolution(size, solution);
}

/** set branching priority */
void setBranchPriorities(
		DspApiEnv * env,       /**< pointer to API object */
		int         size,      /**< number of priorities */
		int *       priorities /**< A smaller number has more priority. */)
{
	getTssModel(env)->setPriorities(size, priorities);
}

/** solve deterministic equivalent model */
void solveDe(DspApiEnv * env)
{
	STO_API_CHECK_MODEL();
	freeSolver(env);

	if (!prepareDecModel(env))
		return;

	env->solver_ = new DeDriver(env->par_, env->model_);
	env->solver_->init();
	env->solver_->run();
}

/** solve dual decomposition */
void solveDd(DspApiEnv * env, MPI_Comm comm)
{
	STO_API_CHECK_MODEL();
	freeSolver(env);

	if (!prepareDecModel(env))
		return;

	env->solver_ = new DdDriver(env->par_, env->model_, comm);
	env->solver_->init();
	env->solver_->run();
}

/** solve Benders decomposition using MPI */
void solveBd(
		DspApiEnv * env,      /**< pointer to API object */
		int         nauxvars, /**< number of auxiliary variables (scenario clusters) */
		MPI_Comm    comm      /**< MPI communicator */)
{
	STO_API_CHECK_MODEL();
	freeSolver(env);

	if (!env->model_->isStochastic())
	{
		printf("Error: Benders decomposition is currently only supported by stochastic models\n");
		return;
	}

	if (comm == MPI_UNDEFINED)
		env->solver_ = new BdDriver(env->par_, new DecTssModel(*getTssModel(env)));
	else
		env->solver_ = new BdDriver(env->par_, new DecTssModel(*getTssModel(env)), comm);
	BdDriver * bd = dynamic_cast<BdDriver*>(env->solver_);

	double * obj_aux  = NULL;
	double * clbd_aux = NULL;
	double * cubd_aux = NULL;

	if (nauxvars > 1)
	{
		/** set auxiliary variables */
		obj_aux  = new double [nauxvars];
		clbd_aux = new double [nauxvars];
		cubd_aux = new double [nauxvars];
		CoinFillN(obj_aux, nauxvars, 1.0);
		CoinFillN(clbd_aux, nauxvars, -COIN_DBL_MAX);
		CoinFillN(cubd_aux, nauxvars, +COIN_DBL_MAX);

		bd->setAuxVarData(nauxvars, obj_aux, clbd_aux, cubd_aux);

		FREE_ARRAY_PTR(obj_aux);
		FREE_ARRAY_PTR(clbd_aux);
		FREE_ARRAY_PTR(cubd_aux);
	}

	/** relax second-stage integrality */
	env->par_->setBoolPtrParam("RELAX_INTEGRALITY", 1, true);

	env->solver_->init();
	env->solver_->run();
}

/** set boolean parameter */
void setBoolParam(DspApiEnv * env, const char * name, bool value)
{
	STO_API_CHECK_ENV();
	string strname(name);
	env->par_->setBoolParam(name, value);
}

/** set integer parameter */
void setIntParam(DspApiEnv * env, const char * name, int value)
{
	STO_API_CHECK_ENV();
	string strname(name);
	env->par_->setIntParam(strname, value);
}

/** set double parameter */
void setDblParam(DspApiEnv * env, const char * name, double value)
{
	STO_API_CHECK_ENV();
	string strname(name);
	env->par_->setDblParam(strname, value);
}

/** set string parameter */
void setStrParam(DspApiEnv * env, const char * name, const char *  value)
{
	STO_API_CHECK_ENV();
	string strname(name);
	string strvalue(value);
	env->par_->setStrParam(strname, strvalue);
}

/** set boolean pointer parameter */
void setBoolPtrParam(DspApiEnv * env, const char * name, int size, bool * value)
{
	STO_API_CHECK_ENV();
	string strname(name);
	env->par_->setBoolPtrParamSize(strname, size);
	for (int i = 0; i < size; ++i)
		env->par_->setBoolPtrParam(strname, i, value[i]);
}

/** set integer pointer parameter */
void setIntPtrParam(DspApiEnv * env, const char * name, int size, int * value)
{
	STO_API_CHECK_ENV();
	string strname(name);
	env->par_->setIntPtrParamSize(strname, size);
	for (int i = 0; i < size; ++i)
		env->par_->setIntPtrParam(strname, i, value[i]);
}

/**
 * Get functions
 */

/** get number of rows */
int getNumRows(DspApiEnv * env, int stage)
{
	STO_API_CHECK_MODEL(-1);
	return getTssModel(env)->getNumRows(stage);
}

/** get number of columns */
int getNumCols(DspApiEnv * env, int stage)
{
	STO_API_CHECK_MODEL(-1);
	return getTssModel(env)->getNumCols(stage);
}

/** get number of scenarios */
int getNumScenarios(DspApiEnv * env)
{
	STO_API_CHECK_MODEL(-1);
	return getTssModel(env)->getNumScenarios();
}

/** get number of columns */
int getTotalNumCols(DspApiEnv * env)
{
	STO_API_CHECK_MODEL(-1);
	return env->model_->getFullModelNumCols();
}

/** get number of subproblems */
int getNumSubproblems(DspApiEnv * env)
{
	STO_API_CHECK_MODEL(-1);
	return env->model_->getNumSubproblems();
}

void getObjCoef(DspApiEnv * env, double * obj)
{
	STO_API_CHECK_SOLVER();
	return env->model_->getObjCoef(obj);
}

/** get total cpu time */
double getCpuTime(DspApiEnv * env)
{
	STO_API_CHECK_SOLVER(0.0);
	return env->solver_->getCpuTime();
}

/** get solution time */
double getWallTime(DspApiEnv * env)
{
	STO_API_CHECK_SOLVER(0.);
	return env->solver_->getWallTime();
}

/** get solution status */
int getStatus(DspApiEnv * env)
{
	STO_API_CHECK_SOLVER(0);
	return env->solver_->getStatus();
}

/** get objective value */
double getPrimalBound(DspApiEnv * env)
{
	STO_API_CHECK_SOLVER(0.0);
	return env->solver_->getPrimalObjectiveValue();
}

/** get objective value */
double getDualBound(DspApiEnv * env)
{
	STO_API_CHECK_SOLVER(0.0);
	return env->solver_->getDualObjectiveValue();
}

/** get solution */
void getSolution(DspApiEnv * env, int num, double * solution)
{
	STO_API_CHECK_SOLVER();
	CoinCopyN(env->solver_->getPrimalSolution(), num, solution);
}

/** get number of iterations */
int getNumIterations(DspApiEnv * env)
{
	STO_API_CHECK_SOLVER(0);
	return env->solver_->getNumIterations();
}

/** get number of nodes */
int getNumNodes(DspApiEnv * env)
{
	STO_API_CHECK_SOLVER(0);
	return env->solver_->getNumNodes();
}

/** get number of infeasible solutions evaluated in DD */
int getDdNumInfeasSolutions(DspApiEnv * env)
{
	STO_API_CHECK_SOLVER(0);
	DdDriver * dd = dynamic_cast<DdDriver*>(env->solver_);
	if (dd)
		return dd->getNumInfeasibleSolutions();
	return 0;
}

/** get number of master problem solved in DD */
int getDdNumMasterSolved(DspApiEnv * env)
{
	STO_API_CHECK_SOLVER(0);
	DdDriver * dd = dynamic_cast<DdDriver*>(env->solver_);
	if (dd)
		return dd->getMasterPtr()->s_statuses_.size();
	return 0;
}

/** get master solution cpu time per iteration in DD */
void getDdMasterCpuTimes(DspApiEnv * env, double * time)
{
	STO_API_CHECK_SOLVER();
	DdDriver * dd = dynamic_cast<DdDriver*>(env->solver_);
	if (dd)
	{
		for (unsigned i = 0; i < dd->getMasterPtr()->s_cputimes_.size(); ++i)
			time[i] = dd->getMasterPtr()->s_cputimes_[i];
	}
}

/** get master solution wall time per iteration in DD */
void getDdMasterWallTimes(DspApiEnv * env, double * time)
{
	STO_API_CHECK_SOLVER();
	DdDriver * dd = dynamic_cast<DdDriver*>(env->solver_);
	if (dd)
	{
		for (unsigned i = 0; i < dd->getMasterPtr()->s_walltimes_.size(); ++i)
			time[i] = dd->getMasterPtr()->s_walltimes_[i];
	}
}

/** get history of master primal objective values */
void getDdMasterPrimalBounds(DspApiEnv * env, double * vals)
{
	STO_API_CHECK_SOLVER();
	DdDriver * dd = dynamic_cast<DdDriver*>(env->solver_);
	if (dd)
	{
		for (unsigned i = 0; i < dd->getMasterPtr()->s_primobjs_.size(); ++i)
			vals[i] = dd->getMasterPtr()->s_primobjs_[i];
	}
}

/** get history of master dual objective values */
void getDdMasterDualBounds(DspApiEnv * env, double * vals)
{
	STO_API_CHECK_SOLVER();
	DdDriver * dd = dynamic_cast<DdDriver*>(env->solver_);
	if (dd)
	{
		for (unsigned i = 0; i < dd->getMasterPtr()->s_dualobjs_.size(); ++i)
			vals[i] = dd->getMasterPtr()->s_dualobjs_[i];
	}
}

/** get number of subproblems solved in DD */
int getDdNumSubproblemSolved(DspApiEnv * env)
{
	STO_API_CHECK_SOLVER(0);
	DdDriver * dd = dynamic_cast<DdDriver*>(env->solver_);
	if (dd)
		return  dd->getWorkerPtr()->s_statuses_.size();
	return 0;
}

/** get master primal solutionshistory of dual multipliers from dual decomposition */
void getDdMasterPrimalSolutions(
		DspApiEnv* env, /**< API environment */
		double**   vals /**< array of multipliers */)
{
	STO_API_CHECK_SOLVER();
	DdDriver * dd = dynamic_cast<DdDriver*>(env->solver_);
	if (dd)
	{
		for (unsigned i = 0; i < dd->getMasterPtr()->s_primsols_.size(); ++i)
			for (int j = 0; j < env->model_->getNumCouplingRows(); ++j)
				vals[i][j] = dd->getMasterPtr()->s_primsols_[i][j];
	}
}

/** get solution time per cpu iteration in DD */
void getDdSubproblemCpuTimes(DspApiEnv * env, double * time)
{
	STO_API_CHECK_SOLVER();
	DdDriver * dd = dynamic_cast<DdDriver*>(env->solver_);
	if (dd)
	{
		for (unsigned i = 0; i < dd->getWorkerPtr()->s_cputimes_.size(); ++i)
			time[i] = dd->getWorkerPtr()->s_cputimes_[i];
	}
}

/** get solution time per wall iteration in DD */
void getDdSubproblemWallTimes(DspApiEnv * env, double * time)
{
	STO_API_CHECK_SOLVER();
	DdDriver * dd = dynamic_cast<DdDriver*>(env->solver_);
	if (dd)
	{
		for (unsigned i = 0; i < dd->getWorkerPtr()->s_walltimes_.size(); ++i)
			time[i] = dd->getWorkerPtr()->s_walltimes_[i];
	}
}

/** get history of subproblem primal objective values */
void getDdSubproblemPrimalBounds(DspApiEnv * env, double * vals)
{
	STO_API_CHECK_SOLVER();
	DdDriver * dd = dynamic_cast<DdDriver*>(env->solver_);
	if (dd)
	{
		for (unsigned i = 0; i < dd->getWorkerPtr()->s_primobjs_.size(); ++i)
			vals[i] = dd->getWorkerPtr()->s_primobjs_[i];
	}
}

/** get history of subproblem dual objective values */
void getDdSubproblemDualBounds(DspApiEnv * env, double * vals)
{
	STO_API_CHECK_SOLVER();
	DdDriver * dd = dynamic_cast<DdDriver*>(env->solver_);
	if (dd)
	{
		for (unsigned i = 0; i < dd->getWorkerPtr()->s_dualobjs_.size(); ++i)
			vals[i] = dd->getWorkerPtr()->s_dualobjs_[i];
	}
}

/**
 * Misc
 */

/** print model */
void printModel(DspApiEnv * env)
{
	STO_API_CHECK_MODEL();
	env->model_->__printData();
}

#undef STO_API_CHECK_MODEL
#undef STO_API_CHECK_ENV

#ifdef __cplusplus
}

#endif


