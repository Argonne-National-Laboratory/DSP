/*
 * TssBdSub.h
 *
 *  Created on: Dec 5, 2014
 *      Author: kibaekkim
 */

#ifndef TSSBDSUB_H_
#define TSSBDSUB_H_

#include <assert.h>

/** Open MP */
#ifdef USE_OMP
#include <omp.h>
#endif

/** Coin-OR */
#include "OsiSolverInterface.hpp"
#include "OsiCuts.hpp"
#include "CoinTime.hpp"

/** Dsp */
#include "Model/TssModel.h"
#include "Utility/StoConfig.h"
#include "Utility/StoMacros.h"
#include "Utility/StoRtnCodes.h"
#include "Utility/StoUtility.h"
#include "SolverInterface/SolverInterface.h"
#include "Solver/StoParam.h"

#define DSP_SOLVER_INTERFACE

class TssBdSub
{
public:

	enum whereFrom
	{
		TssBd = 0,
		TssDd
	};

	enum whichCut
	{
		BothCuts = 0,
		FeasCut,
		OptCut
	};

	/** default constructor */
	TssBdSub(StoParam * par);

	/** default destructor */
	virtual ~TssBdSub();

	/** load problem */
	STO_RTN_CODE loadProblem(
			TssModel * model,
			int naugs    = 0,
			int * augs   = NULL,
			int nAuxvars = 1);

#if 0
	/** load problem */
	STO_RTN_CODE loadProblem(
			int nSubs,
			const double * weight,
			const CoinPackedMatrix ** mat_mp,
			const OsiSolverInterface ** cglp,
			int nAuxvars = 1);
#endif

	/** generate Benders cut */
	void generateCuts(
			int            ncols,
			const double * x,
			OsiCuts *      cs,
			int            where = TssBd,
			int            which = BothCuts);

	/** solve single scenario recourse */
	void solveSingleRecourse(
			int            scenario, /**< scenario index */
			const double * x,        /**< first-stage variables */
			double &       objval,   /**< second-stage objective value */
			double *       solution  /**< second-stage solution */);

	/** solve recourse (without cut generation) */
	void solveRecourse(
			const double * x,         /**< first-stage variables */
			double *       objval,    /**< second-stage objective values */
			double **      solution,  /**< second-stage solutions */
			int            ncores = 1 /**< number of cores used to run in parallel */);

private:

	/** solve one subproblem. this is a body of loop in gutsOfGenerateCuts */
	static void solveOneFeasSubproblem(
			TssBdSub *     cgl,
			int            s,      /**< scenario index */
			const double * x,      /**< first-stage solution */
			double **      Tx,     /**< Tx */
			double **      cutval, /**< Benders cut body */
			double *       cutrhs  /**< Benders cut RHS */);

	/** solve one subproblem. this is a body of loop in gutsOfGenerateCuts */
	static void solveOneOptSubproblem(
			TssBdSub *     cgl,
			int            s,      /**< scenario index */
			const double * x,      /**< first-stage solution */
			double **      Tx,     /**< Tx */
			double **      cutval, /**< Benders cut body */
			double *       cutrhs  /**< Benders cut RHS */);

	/** solve one subproblem. this is a body of loop in gutsOfGenerateCuts */
	static void solveOneSubproblem(
			TssBdSub *     cgl,
			int            s,                /**< scenario index */
			const double * x,                /**< first-stage solution */
			double **      Tx,               /**< Tx */
			double **      cutval,           /**< Benders cut body */
			double *       cutrhs,           /**< Benders cut RHS */
			int            enableOptCuts = 1 /**< whether to generate optimality cuts or not */);

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

	/** construct cuts (e.g. aggregation) */
	int constructCuts(
			int            ncols,         /**< [in] number of columns in the master */
			double **      values,        /**< [in] dense vector of cut coefficients in the size of nSubs_ */
			double *       rhs,           /**< [in] array of cut rhs in the size of nSubs_ */
			const double * x,             /**< [in] current solution to calculate effectiveness */
			OsiCuts *      cuts           /**< [out] a set of cuts generated */);

	/** solve feasibility problem */
	static int solveFeasProblem(
			OsiSolverInterface * si, /**< [in] subproblem solver interface */
			int & nAddedCols         /**< [out] number of columns added */);

	/** change feasibility problem to original problem */
	static int chgToOrgProblem(
			OsiSolverInterface * si, /**< [in] subproblem solver interface */
			const double * obj,      /**< [in] original objective function */
			int & nAddedCols         /**< [out] number of columns added */);

public:
	StoParam * par_; /**< parameters */

	int nSubs_;                    /**< number of subproblems (CGLPs) */
	int nAuxvars_;                 /**< number of auxiliary variables for Benders decomposition */
	bool * excludedScenarios_;     /**< mark whether scenario is excluded and augmented to the master */
	const double * weight_;        /**< weight to cut for each subproblem */
	CoinPackedMatrix   ** mat_mp_; /**< array of matrix corresponding to master problem part */
	SolverInterface ** recourse_;  /**< array of solver interfaces */
	OsiSolverInterface ** cglp_;   /**< array of Cut Generation LP */
	CoinWarmStart ** warm_start_;  /**< warm start information for each subproblem */

	double cut_generation_wall_time_;    /**< wall time in cut generation */
	int nFeasSolutions_; /**< number of feasible solutions visited */

	vector<int> scenarios_; /**< scenario indices for cut generation */

	STO_RTN_CODE * status_; /** subproblem solution status */
};

#endif /* TSSBDSUB_H_ */
