/*
 * DspCInterface.h
 *
 *  Created on: Aug 30, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_DSPCINTERFACE_H_
#define SRC_DSPCINTERFACE_H_

/** standard libraries */
#ifdef DSP_HAS_MPI
#include "mpi.h"
#endif
/** coin libraries */
#include "CoinHelperFunctions.hpp"
/** Dsp libraries */
#include "DspApiEnv.h"

/** forward definition */
class TssModel;

/* using __cplusplus */
#ifdef __cplusplus
extern "C" {
#endif

/** create API environment */
DspApiEnv * createEnv(void);

/** free API environment */
void freeEnv(DspApiEnv * env);

/** free new model */
void freeModel(DspApiEnv * env);

/** free solver */
void freeSolver(DspApiEnv * env);

/** If current model is stochastic, return the model as a TssModel object. If no model exists, create one. */
TssModel * getTssModel(DspApiEnv * env);

/** get model pointer */
DecModel * getModelPtr(DspApiEnv * env);

/**
 * The following functions are specialized for stochastic programming.
 */

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
void readSmps(DspApiEnv * env, const char * smps);

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

/**
 * The following two functions provide a way to store a structured problem. The first function
 * stores a deterministic equivalent problem (a.k.a. extensive form), and then the second function
 * specifies number of subproblems, coupling constraints, and subproblem columns. However, this does
 * not consider distributed memory.
 */

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
		const double *       rubd    /**< row upper bounds */);

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
		double *    couplingRhs     /**< right-hand sides of each coupling constraint */);

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

/** solve decomposition */
void solveDw(DspApiEnv * env);

#ifdef DSP_HAS_MPI
/** solve decomposition */
void solveDwMpi(DspApiEnv * env, MPI_Comm comm);
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
void setStrParam(DspApiEnv * env, const char * name, const char *  value);

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

/** get number of scenarios */
int getNumScenarios(DspApiEnv * env);

/** get number of coupling rows */
int getNumCouplingRows(DspApiEnv * env);

/** get total number of rows */
int getTotalNumRows(DspApiEnv * env);

/** get total number of columns */
int getTotalNumCols(DspApiEnv * env);

/** get number of subproblems */
int getNumSubproblems(DspApiEnv * env);

void getObjCoef(DspApiEnv * env, double * obj);

/** get total cpu time */
double getCpuTime(DspApiEnv * env);

/** get solution time */
double getWallTime(DspApiEnv * env);

/** get solution status */
int getStatus(DspApiEnv * env);

/** get objective value */
double getPrimalBound(DspApiEnv * env);

/** get objective value */
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
 * Misc
 */

void writeMps(DspApiEnv * env, const char * name);

/** print model */
void printModel(DspApiEnv * env);

#ifdef __cplusplus
}
#endif

#endif /* SRC_DSPCINTERFACE_H_ */
