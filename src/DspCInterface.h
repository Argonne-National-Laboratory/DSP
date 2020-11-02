/*
 * DspCInterface.h
 *
 *  Created on: Oct 24, 2014
 *      Author: kibaekkim
 */

#ifndef DSPCINTERFACE_H_
#define DSPCINTERFACE_H_

#include <cstdlib>
#include <cstdio>

#include "CoinTypes.hpp"
#include "Utility/DspMpi.h"

/** DSP */
#include "DspApiEnv.h"

class TssModel;

#ifdef __cplusplus
extern "C" {
#endif

/** create API environment */
DspApiEnv * createEnv(void);

/** free API environment */
void freeEnv(DspApiEnv * env);

/** free model */
void freeModel(DspApiEnv * env);

/** free solver */
void freeSolver(DspApiEnv * env);

/** If current model is stochastic, return the model as a TssModel object. If no model exists, create one. */
TssModel * getTssModel(DspApiEnv * env);

/** get model pointer */
DecModel * getModelPtr(DspApiEnv * env);

/** set number of scenarios */
void setNumberOfScenarios(
		DspApiEnv * env,  /**< pointer to API object */
		int         nscen /**< number of scenarios */);

/** set dimensions */
void setDimensions(
		DspApiEnv * env,    /**< pointer to API object */
		const int   ncols1, /**< number of first-stage columns */
		const int   nrows1, /**< number of first-stage rows */
		const int   ncols2, /**< number of second-stage columns */
		const int   nrows2  /**< number of second-stage rows */);

/** read smps files */
int readSmps(DspApiEnv * env, const char * smps);

/** read dro files */
int readDro(DspApiEnv * env, const char * dro);

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
 		const double *       rubd   /**< row upper bounds */);

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
		const double *       rubd   /**< row upper bounds */);

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
 		const double *       rubd   /**< row upper bounds */);

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
		const double *       rubd   /**< row upper bounds */);		

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
		const double *       rubd    /**< row upper bounds */);

/** update block structure information */
void updateBlocks(
		DspApiEnv * env /**< pointer to API object */);

/** set initial solutions
 * This function can be called multiple times for multiple initial solutions. */
void setSolution(
		DspApiEnv * env,     /**< pointer to API object */
		int         size,    /**< size of solution array */
		double *    solution /**< solution to set */);

/** solve deterministic equivalent model */
void solveDe(DspApiEnv * env /**< pointer to API object */);

/** solve dual decomposition */
void solveDd(DspApiEnv * env /**< pointer to API object */);

/** solve Dantzig-Wolfe decomposition with branch-and-bound */
void solveDw(DspApiEnv * env /**< pointer to API object */);

/** solve serial Benders decomposition */
void solveBd(DspApiEnv * env /**< pointer to API object */);

#ifdef DSP_HAS_MPI
/** solve parallel dual decomposition */
void solveDdMpi(
		DspApiEnv * env, /**< pointer to API object */
		MPI_Comm    comm /**< MPI communicator */);

/** solve parallel Dantzig-Wolfe decomposition with branch-and-bound */
void solveDwMpi(
		DspApiEnv * env, /**< pointer to API object */
		MPI_Comm    comm /**< MPI communicator */);

/** solve parallel Benders decomposition */
void solveBdMpi(
		DspApiEnv * env, /**< pointer to API object */
		MPI_Comm    comm /**< MPI communicator */);
#endif

/** read parameter file */
void readParamFile(DspApiEnv * env, const char * param_file);

/** set boolean parameter */
void setBoolParam(DspApiEnv * env, const char * name, bool value);

/** set integer parameter */
void setIntParam(DspApiEnv * env, const char * name, int value);

/** set double parameter */
void setDblParam(DspApiEnv * env, const char * name, double value);

/** set string parameter */
void setStrParam(DspApiEnv * env, const char * name, const char * value);

/** set boolean pointer parameter */
void setBoolPtrParam(DspApiEnv * env, const char * name, int size, bool * value);

/** set integer pointer parameter */
void setIntPtrParam(DspApiEnv * env, const char * name, int size, int * value);

/**
 * Get functions
 */

/** get number of rows */
int getNumRows(DspApiEnv * env, int stage);

/** get number of columns */
int getNumCols(DspApiEnv * env, int stage);

/** get number of integer variables */
int getNumIntegers(DspApiEnv * env, int stage);

/** get Total number of row */
int getTotalNumRows(DspApiEnv * env);

/** get number of columns */
int getTotalNumCols(DspApiEnv * env);

/** get number of scenarios */
int getNumScenarios(DspApiEnv * env);

/** get number of coupling rows */
int getNumCouplingRows(DspApiEnv * env);

/** get number of subproblems */
int getNumSubproblems(DspApiEnv * env);

/** get wall cpu time */
double getCpuTime(DspApiEnv * env);

/** get solution wall time */
double getWallTime(DspApiEnv * env);

/** get solution status */
int getStatus(DspApiEnv * env);

/** get primal bound */
double getPrimalBound(DspApiEnv * env);

/** get dual bound */
double getDualBound(DspApiEnv * env);

/** get solution */
void getPrimalSolution(DspApiEnv * env, int num, double * solution);

/** get dual solution */
void getDualSolution(DspApiEnv * env, int num, double * solution);

/** get number of iterations */
int getNumIterations(DspApiEnv * env);

/** get number of nodes */
int getNumNodes(DspApiEnv * env);

/**
 * Misc.
 */

/** write extensive form MPS */
void writeMps(DspApiEnv * env, const char * name);

/** print model */
void printModel(DspApiEnv *env);

#ifdef __cplusplus
}
#endif

#endif /* DSPCINTERFACE_H_ */
