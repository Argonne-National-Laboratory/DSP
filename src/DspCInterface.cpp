/*
 * DspCInterface.cpp
 *
 *  Created on: Oct 27, 2014
 *      Author: kibaekkim
 */

// #define DSP_DEBUG

#include <cstdlib>
#include <cstdio>

#include "DspApiEnv.h"
#include "DspCInterface.h"
#include "Utility/DspMacros.h"
#include "Model/DecTssModel.h"
#include "Model/DecBlkModel.h"

#include "Solver/Deterministic/DeDriver.h"
#include "Solver/DualDecomp/DdDriverSerial.h"
#include "Solver/DantzigWolfe/DwSolverSerial.h"

#ifdef DSP_HAS_SCIP
#include "Solver/Benders/BdDriverSerial.h"
#endif /* DSP_HAS_SCIP */


#ifdef DSP_HAS_MPI

#ifdef DSP_HAS_SCIP
#include "Solver/Benders/BdDriverMpi.h"
#endif /* DSP_HAS_SCIP */

#include "Solver/DualDecomp/DdDriverMpi.h"
#include "Solver/DantzigWolfe/DwSolverMpi.h"

#endif /* DSP_HAS_MPI */



/* using __cplusplus */
#ifdef __cplusplus
extern "C" {
#endif

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
	DSP_API_CHECK_ENV();
	FREE_PTR(env->model_);
}

/** free solver */
void freeSolver(DspApiEnv * env)
{
	DSP_API_CHECK_ENV();
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

/** get model pointer */
DecModel * getModelPtr(DspApiEnv * env)
{
	DSP_API_CHECK_ENV(NULL);
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
int readSmps(DspApiEnv * env, const char * smps)
{
	return getTssModel(env)->readSmps(smps);
}

/** read dro files */
int readDro(DspApiEnv * env, const char * dro)
{
	return getTssModel(env)->readDro(dro);
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

/** load first-stage problem with quadratic objective */
void loadQuadraticFirstStage(
		DspApiEnv *          env,   /**< pointer to API object */
		const CoinBigIndex * start, /**< start index for each row */
		const int *          index, /**< column indices */
		const double *       value, /**< constraint elements */
		const double *       clbd,  /**< column lower bounds */
		const double *       cubd,  /**< column upper bounds */
		const char *         ctype, /**< column types */
		const double *       obj,   /**< objective coefficients */
		const int * 		 qrowindex,/**< start index for first-stage quadratic objective */
		const int *			 qcolindex,/**< quadratic objective column indices */
		const double *		 qvalue,/**< quadratic objective elements */
		const CoinBigIndex	 qnum,	/**< number of quadratic elements */
		const double *       rlbd,  /**< row lower bounds */
		const double *       rubd   /**< row upper bounds */)
{
	getTssModel(env)->loadFirstStage(start, index, value, clbd, cubd, ctype, obj, qrowindex, qcolindex, qvalue, qnum, rlbd, rubd);
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

/** load second-stage problem */
void loadQuadraticSecondStage(
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
		const int * 		 qrowindex,/**< start index for first-stage quadratic objective */
		const int *			 qcolindex,/**< quadratic objective column indices */
		const double *		 qvalue,/**< quadratic objective elements */
		const CoinBigIndex	 qnum,	/**< number of quadratic elements */
		const double *       rlbd,  /**< row lower bounds */
		const double *       rubd   /**< row upper bounds */)
{
	getTssModel(env)->loadSecondStage(s, prob, start, index, value, clbd, cubd, ctype, obj, qrowindex, qcolindex, qvalue, qnum, rlbd, rubd);
}

/**
 * The following function allows to read the model by blocks. This
 * is particularly useful for distributed memory computing.
 */

/** load block problems */
void loadBlockProblem(
		DspApiEnv *          env,    /**< pointer to API object */
		int                  id,     /**< block ID */
		int                  ncols,  /**< number of columns */
		int                  nrows,  /**< number of rows */
		int                  numels, /**< number of elements in the matrix */
		const CoinBigIndex * start,  /**< start index for each row */
		const int *          index,  /**< column indices */
		const double *       value,  /**< constraint elements */
		const double *       clbd,   /**< column lower bounds */
		const double *       cubd,   /**< column upper bounds */
		const char *         ctype,  /**< column types */
		const double *       obj,    /**< objective coefficients */
		const double *       rlbd,   /**< row lower bounds */
		const double *       rubd    /**< row upper bounds */) {
	DSPdebugMessage("\n##############################################################\n"
			          "## The beginning of loadBlockProblem                        ##\n"
			          "##############################################################\n");
	DSPdebugMessage("id = %d\n", id);
	DSPdebugMessage("ncols = %d\n", ncols);
	DSPdebugMessage("nrows = %d\n", nrows);
	DSPdebugMessage("numels = %d\n", numels);
	DSPdebugMessage("start:\n");
	DSPdebug(DspMessage::printArray(nrows+1, start));
	DSPdebugMessage("index:\n");
	DSPdebug(DspMessage::printArray(numels, index));
	DSPdebugMessage("value:\n");
	DSPdebug(DspMessage::printArray(numels, value));
	DSPdebugMessage("clbd:\n");
	DSPdebug(DspMessage::printArray(ncols, clbd));
	DSPdebugMessage("cubd:\n");
	DSPdebug(DspMessage::printArray(ncols, cubd));
	DSPdebugMessage("obj:\n");
	DSPdebug(DspMessage::printArray(ncols, obj));
	DSPdebugMessage("rlbd:\n");
	DSPdebug(DspMessage::printArray(nrows, rlbd));
	DSPdebugMessage("rubd:\n");
	DSPdebug(DspMessage::printArray(nrows, rubd));
	DSPdebugMessage("\n##############################################################\n"
			          "## The end of loadBlockProblem                              ##\n"
			          "##############################################################\n\n\n\n");
	if (env->model_ == NULL)
		env->model_ = new DecBlkModel;
	BlkModel* blk = dynamic_cast<DecBlkModel*>(env->model_)->blkPtr();
	if (blk == NULL)
		fprintf(stderr, "Block decomposition model is not loaded.\n");
	else {
		blk->addBlock(id,
				new DetBlock(start, index, value, numels,
						ncols, nrows, clbd, cubd, ctype, obj, rlbd, rubd));
	}
}

/** update block structure information */
void updateBlocks(
		DspApiEnv * env /**< pointer to API object */) {
	if (env->model_ == NULL)
		fprintf(stderr, "Block decomposition model is not loaded.\n");
	BlkModel* blk = dynamic_cast<DecBlkModel*>(env->model_)->blkPtr();
	if (blk == NULL)
		fprintf(stderr, "Block decomposition model is not loaded.\n");
	else
		blk->updateBlocks();
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

/** solve deterministic equivalent model */
void solveDe(DspApiEnv * env)
{
	BGN_TRY_CATCH

	DSP_API_CHECK_MODEL();
	freeSolver(env);

	env->solver_ = new DeDriver(env->model_, env->par_, env->message_);
	DSP_RTN_CHECK_THROW(env->solver_->init());
	DSP_RTN_CHECK_THROW(dynamic_cast<DeDriver*>(env->solver_)->run());
	DSP_RTN_CHECK_THROW(env->solver_->finalize());

	END_TRY_CATCH(;)
}

/** solve dual decomposition */
void solveDd(DspApiEnv * env)
{
	BGN_TRY_CATCH

	DSP_API_CHECK_MODEL();
	freeSolver(env);

	env->solver_ = new DdDriverSerial(env->model_, env->par_, env->message_);
	
	DSP_RTN_CHECK_THROW(env->solver_->init());
	DSP_RTN_CHECK_THROW(dynamic_cast<DdDriverSerial*>(env->solver_)->run());
	DSP_RTN_CHECK_THROW(env->solver_->finalize());

	END_TRY_CATCH(;)
}

/** solve Dantzig-Wolfe decomposition (with branch-and-bound) */
void solveDw(DspApiEnv * env)
{
	BGN_TRY_CATCH

	DSP_API_CHECK_MODEL();
	freeSolver(env);

	env->solver_ = new DwSolverSerial(env->model_, env->par_, env->message_);
	DSP_RTN_CHECK_THROW(env->solver_->init());
	DSP_RTN_CHECK_THROW(dynamic_cast<DwSolverSerial*>(env->solver_)->run());
	DSP_RTN_CHECK_THROW(env->solver_->finalize());

	END_TRY_CATCH(;)
}

/** solve serial Benders decomposition */
void solveBd(DspApiEnv * env)
{
	BGN_TRY_CATCH
#ifdef DSP_HAS_SCIP
	DSP_API_CHECK_MODEL();
	freeSolver(env);

	if (!env->model_->isStochastic())
	{
		printf("Error: Benders decomposition is currently only supported by stochastic models\n");
		return;
	}

	DecTssModel* dec = new DecTssModel(*getTssModel(env));
	BdDriverSerial * bd = new BdDriverSerial(dec, env->par_, env->message_);
	env->solver_ = bd;
	DSPdebugMessage("Created a serial Benders object\n");

    int nauxvars = env->par_->getIntParam("BD/NUM_CUTS_PER_ITER");
	if (nauxvars > 0)
	{
		double * obj_aux  = NULL;
		double * clbd_aux = NULL;
		double * cubd_aux = NULL;
		/** set auxiliary variables */
		obj_aux  = new double [nauxvars];
		clbd_aux = new double [nauxvars];
		cubd_aux = new double [nauxvars];
		CoinFillN(obj_aux, nauxvars, 1.0);
		CoinFillN(clbd_aux, nauxvars, -COIN_DBL_MAX);
		CoinFillN(cubd_aux, nauxvars, +COIN_DBL_MAX);
		/** set auxiliary variables */
		bd->setAuxVarData(nauxvars, obj_aux, clbd_aux, cubd_aux);
        /** free memory */
		FREE_ARRAY_PTR(obj_aux);
		FREE_ARRAY_PTR(clbd_aux);
		FREE_ARRAY_PTR(cubd_aux);
	}
	DSPdebugMessage("Set auxiliary variable data\n");

	DSP_RTN_CHECK_THROW(env->solver_->init());
	DSP_RTN_CHECK_THROW(dynamic_cast<BdDriverSerial*>(env->solver_)->run());
	DSP_RTN_CHECK_THROW(env->solver_->finalize());
#else
	printf("Benders decomposition has been disabled because SCIP was not available.\n");
#endif
	END_TRY_CATCH(;)
}

#ifdef DSP_HAS_MPI
/** solve parallel dual decomposition */
void solveDdMpi(DspApiEnv * env, MPI_Comm comm)
{
	BGN_TRY_CATCH

	int comm_size;
	MPI_Comm_size(comm, &comm_size);

	if (comm_size == 1) {
		solveDd(env);
		return;
	}

	DSP_API_CHECK_MODEL();
	freeSolver(env);

	env->solver_ = new DdDriverMpi(env->model_, env->par_, env->message_, comm);
	DSP_RTN_CHECK_THROW(env->solver_->init());
	DSP_RTN_CHECK_THROW(dynamic_cast<DdDriverMpi*>(env->solver_)->run());
	DSP_RTN_CHECK_THROW(env->solver_->finalize());

	END_TRY_CATCH(;)
}

/** solve parallel Dantzig-Wolfe decomposition with branch-and-bound */
void solveDwMpi(DspApiEnv * env, MPI_Comm comm)
{
	BGN_TRY_CATCH

	int comm_size;
	MPI_Comm_size(comm, &comm_size);

	if (comm_size == 1) {
		solveDw(env);
		return;
	}

	DSP_API_CHECK_MODEL();
	freeSolver(env);

	env->solver_ = new DwSolverMpi(env->model_, env->par_, env->message_, comm);
	DSP_RTN_CHECK_THROW(env->solver_->init());
	DSP_RTN_CHECK_THROW(dynamic_cast<DwSolverMpi*>(env->solver_)->run());
	DSP_RTN_CHECK_THROW(env->solver_->finalize());

	END_TRY_CATCH(;)
}

/** solve parallel Benders decomposition */
void solveBdMpi(
		DspApiEnv * env, MPI_Comm comm)
{
	BGN_TRY_CATCH
	
	int comm_size;
	MPI_Comm_size(comm, &comm_size);

	if (comm_size == 1) {
		solveBd(env);
		return;
	}

#ifdef DSP_HAS_SCIP
	DSP_API_CHECK_MODEL();
	freeSolver(env);

	if (!env->model_->isStochastic())
	{
		printf("Error: Benders decomposition is currently only supported by stochastic models\n");
		return;
	}

	//DSPdebugMessage("Creating a MPI Benders object (comm %d)\n", comm);
	//BdDriverMpi * bd = new BdDriverMpi(new DecTssModel(*getTssModel(env)), env->par_, env->message_, comm); //this line cause segment fault
	BdDriverMpi * bd = new BdDriverMpi(env->model_, env->par_, env->message_, comm);
	env->solver_ = bd;

	double * obj_aux  = NULL;
	double * clbd_aux = NULL;
	double * cubd_aux = NULL;

    int nauxvars = env->par_->getIntParam("BD/NUM_CUTS_PER_ITER");
	if (nauxvars > 0)
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

	DSP_RTN_CHECK_THROW(env->solver_->init());
	DSP_RTN_CHECK_THROW(dynamic_cast<BdDriverMpi*>(env->solver_)->run());
	DSP_RTN_CHECK_THROW(env->solver_->finalize());
#else
	printf("Benders decomposition has been disabled because SCIP was not available.\n");
#endif

	END_TRY_CATCH(;)
}
#endif

/** read parameter file */
void readParamFile(DspApiEnv * env, const char * param_file)
{
	DSP_API_CHECK_ENV();
	env->par_->readParamFile(param_file);
}

/** set boolean parameter */
void setBoolParam(DspApiEnv * env, const char * name, bool value)
{
	DSP_API_CHECK_ENV();
	string strname(name);
	env->par_->setBoolParam(name, value);
}

/** set integer parameter */
void setIntParam(DspApiEnv * env, const char * name, int value)
{
	DSP_API_CHECK_ENV();
	string strname(name);
	env->par_->setIntParam(strname, value);
}

/** set double parameter */
void setDblParam(DspApiEnv * env, const char * name, double value)
{
	DSP_API_CHECK_ENV();
	string strname(name);
	env->par_->setDblParam(strname, value);
}

/** set string parameter */
void setStrParam(DspApiEnv * env, const char * name, const char *  value)
{
	DSP_API_CHECK_ENV();
	string strname(name);
	string strvalue(value);
	env->par_->setStrParam(strname, strvalue);
}

/** set boolean pointer parameter */
void setBoolPtrParam(DspApiEnv * env, const char * name, int size, bool * value)
{
	DSP_API_CHECK_ENV();
	string strname(name);
	env->par_->setBoolPtrParamSize(strname, size);
	for (int i = 0; i < size; ++i)
		env->par_->setBoolPtrParam(strname, i, value[i]);
}

/** set integer pointer parameter */
void setIntPtrParam(DspApiEnv * env, const char * name, int size, int * value)
{
	DSP_API_CHECK_ENV();
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
	DSP_API_CHECK_MODEL(-1);
	return getTssModel(env)->getNumRows(stage);
}

/** get number of columns */
int getNumCols(DspApiEnv * env, int stage)
{
	DSP_API_CHECK_MODEL(-1);
	return getTssModel(env)->getNumCols(stage);
}

/** get number of integer variables */
int getNumIntegers(DspApiEnv * env, int stage)
{
	DSP_API_CHECK_MODEL(-1);
	return getTssModel(env)->getNumIntegers(stage);
}

/** get number of scenarios */
int getNumScenarios(DspApiEnv * env)
{
	DSP_API_CHECK_MODEL(-1);
	return getTssModel(env)->getNumScenarios();
}

/** get number of coupling rows */
int getNumCouplingRows(DspApiEnv * env)
{
	DSP_API_CHECK_MODEL(-1);
	return getModelPtr(env)->getNumCouplingRows();
}

/** get total number of rows */
int getTotalNumRows(DspApiEnv * env)
{
	DSP_API_CHECK_MODEL(-1);
	return getModelPtr(env)->getFullModelNumRows();
}

/** get total number of columns */
int getTotalNumCols(DspApiEnv * env)
{
	DSP_API_CHECK_MODEL(-1);
	return getModelPtr(env)->getFullModelNumCols();
}

/** get number of subproblems */
int getNumSubproblems(DspApiEnv * env)
{
	DSP_API_CHECK_MODEL(-1);
	return getModelPtr(env)->getNumSubproblems();
}

/** get total cpu time */
double getCpuTime(DspApiEnv * env)
{
	DSP_API_CHECK_SOLVER(0.0);
	return env->solver_->getCpuTime();
}

/** get solution time */
double getWallTime(DspApiEnv * env)
{
	DSP_API_CHECK_SOLVER(0.);
	return env->solver_->getWallTime();
}

/** get solution status */
int getStatus(DspApiEnv * env)
{
	DSP_API_CHECK_SOLVER(0);
	return env->solver_->getStatus();
}

/** get objective value */
double getPrimalBound(DspApiEnv * env)
{
	DSP_API_CHECK_SOLVER(0.0);
	return env->solver_->getBestPrimalObjective();
}

/** get objective value */
double getDualBound(DspApiEnv * env)
{
	DSP_API_CHECK_SOLVER(0.0);
	return env->solver_->getBestDualObjective();
}

/** get solution */
void getPrimalSolution(DspApiEnv * env, int num, double * solution)
{
	DSP_API_CHECK_SOLVER();
	CoinCopyN(env->solver_->getBestPrimalSolution(), num, solution);
}

/** get dual solution */
void getDualSolution(DspApiEnv * env, int num, double * solution)
{
	DSP_API_CHECK_SOLVER();
	CoinCopyN(env->solver_->getBestDualSolution(), num, solution);
}

/** get number of iterations */
int getNumIterations(DspApiEnv * env)
{
	DSP_API_CHECK_SOLVER(0);
	return env->solver_->getNumIterations();
}

/** get number of nodes */
int getNumNodes(DspApiEnv * env)
{
	DSP_API_CHECK_SOLVER(0);
	return env->solver_->getNumNodes();
}

/**
 * Misc
 */

void writeMps(DspApiEnv * env, const char * name)
{
	DSP_API_CHECK_MODEL();
	freeSolver(env);

	env->solver_ = new DeDriver(env->model_, env->par_, env->message_);
	env->solver_->init();
	dynamic_cast<DeDriver*>(env->solver_)->writeExtMps(name);
	env->solver_->finalize();
}

/** print model */
void printModel(DspApiEnv * env)
{
	DSP_API_CHECK_MODEL();
	env->model_->__printData();
}

#undef DSP_API_CHECK_MODEL
#undef DSP_API_CHECK_ENV

#ifdef __cplusplus
}
#endif


