/*
 * StoCInterface.h
 *
 *  Created on: Oct 24, 2014
 *      Author: kibaekkim
 */

#ifndef STOCINTERFACE_H_
#define STOCINTERFACE_H_

#include <cstdlib>
#include <cstdio>
#include "mpi.h"
#include "CoinTypes.hpp"
#include "StoApiEnv.h"

class TssModel;

#ifdef __cplusplus
extern "C" {
#endif

/** create API environment */
StoApiEnv * createEnv(void);

/** free API environment */
void freeEnv(StoApiEnv * env);

/** free model */
void freeTssModel(StoApiEnv * env);

/** free solver */
void freeTssSolver(StoApiEnv * env);

/** get model pointer */
TssModel * getModelPtr(StoApiEnv * env);

/** set number of scenarios */
void setNumberOfScenarios(
		StoApiEnv * env,  /**< pointer to API object */
		int         nscen /**< number of scenarios */);

/** set dimensions */
void setDimensions(
		StoApiEnv * env,    /**< pointer to API object */
		const int   ncols1, /**< number of first-stage columns */
		const int   nrows1, /**< number of first-stage rows */
		const int   ncols2, /**< number of second-stage columns */
		const int   nrows2  /**< number of second-stage rows */);

/** read smps files */
void readSmps(StoApiEnv * env, const char * smps);

/** write the extensive form in Mps file */
void writeMps(StoApiEnv * env, const char * mps);

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
		const double *       rubd   /**< row upper bounds */);

/** load first-stage problem */
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
		const double *       rubd   /**< row upper bounds */);

/** set branching priorities */
void branchPriorities(
		StoApiEnv * env,       /**< pointer to API object */
		int *       priorities /**< A smaller number has more priority.
		                         The size of the array is at least the number of integer variables. */);

#if 0
/** add branching object */
void addBranchingObject(
		StoApiEnv * env,     /**< pointer to API object */
		int         nzcnt,   /**< number of nonzero elements */
		int *       indices, /**< indices */
		double *    values,  /**< hyper-plane coefficients */
		int         priority /**< branching priority */);
#endif

/** evaluate solution */
void evaluateSolution(StoApiEnv * env, double * solution);

/** solve deterministic equivalent model */
void solveDe(StoApiEnv * env);

/** solve dual decomposition */
void solveDd(StoApiEnv * env, MPI_Comm comm);

/** solve Benders decomposition */
void solveBd(
		StoApiEnv * env,        /**< pointer to API object */
		int         nauxvars = 1 /**< number of auxiliary variables (scenario clusters) */);

/** set log level */
void setLogLevel(StoApiEnv * env, int level);

/** set number of cores */
void setNumCores(StoApiEnv * env, int num);

/** set node limit */
void setNodeLimit(StoApiEnv * env, int num);

/** set iteration limit */
void setIterLimit(StoApiEnv * env, int num);

/** set wallclock limit */
void setWallLimit(StoApiEnv * env, double lim);

/** set integrality relaxation */
void setIntRelax(StoApiEnv * env, int stage);

/** set Benders augmented scenarios */
void setBdAugScenarios(StoApiEnv * env, int size, int * scenarios);

/** set Benders aggressiveness */
void setBendersAggressive(StoApiEnv * env, int aggressive);

/** set a set of scenarios for the current process */
void setDdProcIdxSet(StoApiEnv * env, int size, int * scenarios);

/** set parameter for adding feasibility cuts */
void setDdAddFeasCuts(StoApiEnv * env, int freq);

/** set parameter for adding optimality cuts */
void setDdAddOptCuts(StoApiEnv * env, int freq);

/** set parameter for evaluating upper bound */
void setDdEvalUb(StoApiEnv * env, int freq);

/** set parameter for logging changes of distance of dual variables */
void setDdDualVarsLog(StoApiEnv * env, int yesNo);

/** set on/off DD recourse cache */
void setDdCacheRecourse(StoApiEnv * env, int yesNo);

/** set Lagrangian master solver */
void setDdMasterSolver(StoApiEnv * env, int type);

/** set DD stopping tolerance */
void setDdStoppingTolerance(StoApiEnv * env, double tol);

/** set number of cuts per iteration added to master */
void setDdMasterNumCutsPerIter(StoApiEnv* env, int num);

/** set maximum number of solutions to evalute */
void setDdMaxPrimsolEval(StoApiEnv* env, int num);

/** set SCIP/display/freq */
void setScipDisplayFreq(StoApiEnv * env, int freq);

/** set SCIP/limits/gap */
void setScipLimitsGap(StoApiEnv * env, double gap);

/** set SCIP/limits/time */
void setScipLimitsTime(StoApiEnv * env, double time);

/**
 * Get functions
 */

/** get number of rows */
int getNumRows(StoApiEnv * env, int stage);

/** get number of columns */
int getNumCols(StoApiEnv * env, int stage);

/** get number of scenarios */
int getNumScenarios(StoApiEnv * env);

/** get solution time */
double getSolutionTime(StoApiEnv * env);

/** get solution status */
int getSolutionStatus(StoApiEnv * env);

void getObjCoef(StoApiEnv * env, double * obj);

/** get objective value */
double getObjValue(StoApiEnv * env);

/** get primal bound */
double getPrimalBound(StoApiEnv * env);

/** get dual bound */
double getDualBound(StoApiEnv * env);

/** get solution */
void getSolution(StoApiEnv * env, int num, double * solution);

/** get number of iterations */
int getNumIterations(StoApiEnv * env);

/** get number of nodes */
int getNumNodes(StoApiEnv * env);

/** get number of infeasible solutions evaluated in DD */
int getDdNumInfeasSolutions(StoApiEnv * env);

/** get size of the iteration time array in DD */
int getDdIterTimeSize(StoApiEnv * env);

/** get iteration time per iteration in DD */
void getDdIterTime(StoApiEnv * env, double * time);

/** get size of the master solution time array in DD */
int getDdMasterTimeSize(StoApiEnv * env);

/** get master solution time per iteration in DD */
void getDdMasterTime(StoApiEnv * env, double * time);

/** get size of the subproblem solution time array in DD */
int getDdSubprobTimeSize(StoApiEnv * env);

/** get subproblem solution time per iteration in DD */
void getDdSubprobTime(StoApiEnv * env, double * time);

/** get number of master objective values */
int getDdNumMasterObjValues(StoApiEnv * env);

/** get number of subproblem objective values */
int getDdNumSubproblemObjValues(StoApiEnv * env);

/** get number of primal bounds */
int getDdNumPrimalBounds(StoApiEnv * env);

/** get number of dual bounds */
int getDdNumDualBounds(StoApiEnv * env);

/** get history of master objective values */
void getDdMasterObjValues(StoApiEnv * env, double * vals);

/** get history of subproblem objective values */
void getDdSubproblemObjValues(StoApiEnv * env, double * vals);

/** get history of primal bounds */
void getDdPrimalBounds(StoApiEnv * env, double * vals);

/** get history of dual bounds */
void getDdDualBounds(StoApiEnv * env, double * vals);

/** get total cpu time */
double getDdCpuTime(StoApiEnv * env);

/** get number of changes of multiplier */
int getDdNumChangesOfMultiplier(StoApiEnv * env);

/** get changes of multiplier */
void getDdChangesOfMultiplier(StoApiEnv * env, double * changes);

/** get DD master time */
double getDdMasterTotalTime(StoApiEnv * env);

/** get DD lower bounding time */
double getDdLbTotalTime(StoApiEnv * env);

/** get DD upper bounding time */
double getDdUbTotalTime(StoApiEnv * env);

/** get DD cut generation time */
double getDdCgTotalTime(StoApiEnv * env);

/**
 * Misc.
 */

/** print model */
void printModel(StoApiEnv *env);

#ifdef __cplusplus
}
#endif

#endif /* STOCINTERFACE_H_ */
