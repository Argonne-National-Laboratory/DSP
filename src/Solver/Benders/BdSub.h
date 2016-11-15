/*
 * BdSub.h
 *
 *  Created on: Feb 18, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_BENDERS_BDSUB_H_
#define SRC_SOLVER_BENDERS_BDSUB_H_

#include <assert.h>

/** DSP */
#include "Utility/DspMacros.h"
#include "Utility/DspRtnCodes.h"
#include "Utility/DspParams.h"
#include "Model/DecModel.h"

/**
 * This class defines a Benders cut generator.
 */
class BdSub {
public:

	/** constructor */
	BdSub(DspParams * par);

	/** destructor */
	virtual ~BdSub();

	/** set subproblem indices */
	DSP_RTN_CODE setSubIndices(int size, int * indices);

	/** load problem */
	DSP_RTN_CODE loadProblem(DecModel * model);

	/** generate Benders cut in raw format (without constructing) */
	int generateCuts(
			int            ncols,  /**< number of first-stage variables */
			const double * x,      /**< first-stage solution */
			double **      cutval, /** dense cut coefficients for each subproblem */
			double *       cutrhs  /** cut rhs for each subproblem */);

public:

	/** get objective value */
	double getObjective(int i) {return objvals_[i];}

	/** get solution */
	const double * getSolution(int i) {return solutions_[i];}

	/** get status for subproblem subindices_[i] */
	DSP_RTN_CODE getStatus(int i) {return status_[i];}

	/** get statuses */
	const int* getStatuses() {return status_;}

	/** get number of subproblems */
	int getNumSubprobs() {return nsubprobs_;}

	/** get subproblem index */
	int getSubprobIndex(int i) {return subindices_[i];}

	/** get subproblem indices */
	int * getSubprobIndices() {return subindices_;}

	/** get number of columns */
	int getNumCols(int i) {return cglp_[i]->getNumCols();}

private:

	/** solve one subproblem. this is a body of loop in gutsOfGenerateCuts */
	static void solveOneSubproblem(
			BdSub *     cgl,
			int            s,                /**< scenario index */
			const double * x,                /**< first-stage solution */
			double **      Tx,               /**< Tx */
			double **      cutval,           /**< Benders cut body */
			double *       cutrhs,           /**< Benders cut RHS */
			int            enableOptCuts = 1 /**< whether to generate optimality cuts or not */);

	/** solve feasibility problem */
	static DSP_RTN_CODE solveFeasProblem(
			OsiSolverInterface * si, /**< [in] subproblem solver interface */
			int & nAddedCols         /**< [out] number of columns added */);

	/** change feasibility problem to original problem */
	static DSP_RTN_CODE chgToOrgProblem(
			OsiSolverInterface * si, /**< [in] subproblem solver interface */
			const double * obj,      /**< [in] original objective function */
			int & nAddedCols         /**< [out] number of columns added */);

	/** calculate cut elements and rhs */
	static int calculateCutElements(
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
			double &       cutrhs              /**< [out] cut rhs */);

private:

	DspParams * par_; /**< parameters */

	int   nsubprobs_;  /**< number of subproblems */
	int * subindices_; /**< subproblem indices indices for cut generation */

	CoinPackedMatrix   ** mat_mp_;     /**< array of matrix corresponding to master problem part */
	OsiSolverInterface ** cglp_;       /**< array of Cut Generation LP */
	CoinWarmStart **      warm_start_; /**< warm start information for each subproblem */
	double *              objvals_;    /**< subproblem objective values */
	double **             solutions_;  /**< subproblem solutions */
	DSP_RTN_CODE *        status_;     /**< subproblem solution status */
};

#endif /* SRC_SOLVER_BENDERS_BDSUB_H_ */
