/*
 * DwSolver.h
 *
 *  Created on: Oct 27, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DANTZIGWOLFE_DWALGO_H_
#define SRC_SOLVER_DANTZIGWOLFE_DWALGO_H_

#include <map>
#include "CoinWarmStartBasis.hpp"
#include "Solver/DecSolver.h"
#include "Solver/DantzigWolfe/DwCol.h"
#include "Solver/DantzigWolfe/DwWorker.h"

/**
 * This creates the master problem data and solves the master problem.
 *
 * The master problem structure is
 *   lb^0 <= \sum_{j=1}^{r^1} H^1 x_j^1 lambda_j^1 + ... + \sum_{j=1}^{r^k} H^k x_j^k lambda_j^k <= ub^0
 *      1 <= \sum_{j=1}^{r^1}           lambda_j^1                                               <= 1
 *    ...
 *      1 <=                                               \sum_{j=1}^{r^k}           lambda_j^k <= 1
 *      0 <=                            lambda_j^1, ...                               lambda_j^k
 *
 * TODO: Branching may be performed on \sum_{j=1}^{r^k} lambda_j^k.
 */

/**
 * TODO: Another way is considering the master problem in the form
 *   lb^0 <= \sum_{j=1}^r (\sum_{i=1}^k H^i x_j^i) lambda_j <= ub^0
 *      1 <= \sum_{j=1}^r                          lambda_j <= 1
 *
 * But not implemented here.
 */
class DwAlgo: public DecSolver {
public:
    /** constructor with worker */
	DwAlgo(DwWorker* worker);

    /** default destructor */
	virtual ~DwAlgo();

    /** initialize */
    virtual DSP_RTN_CODE init() {return DSP_RTN_OK;}

    /** solve */
    virtual DSP_RTN_CODE solve();

    /** finalize */
    virtual DSP_RTN_CODE finalize() {return DSP_RTN_OK;}

protected:

	/** Generate initial columns */
	virtual DSP_RTN_CODE initialColumns();

    /** This creates a master problem. */
    virtual DSP_RTN_CODE createProblem() {return DSP_RTN_OK;}

    /** solve phase 1 */
    virtual DSP_RTN_CODE solvePhase1();

    /** solve phase 2 */
    virtual DSP_RTN_CODE solvePhase2();

    /** guts of solve */
    virtual DSP_RTN_CODE gutsOfSolve();

    /** calculate primal objective value */
    virtual DSP_RTN_CODE calculatePrimalObjective();

    /** restore columns: adding all the columns back */
    virtual DSP_RTN_CODE restoreCols();

    /** reduce columns (e.g., reduced cost fixing) */
    virtual DSP_RTN_CODE reduceCols();

    /** generate columns */
    virtual DSP_RTN_CODE generateCols(
    		const double* price, /**< [in] price */
			double*& piA,        /**< [out] pi^T A */
			double& lb,          /**< [out] lower bound (only for phase 2) */
			int& ncols           /**< [out] number of columns generated */);

    /** calculate piA */
    virtual DSP_RTN_CODE calculatePiA(
    		const double* price, /**< [in] price */
			double*& piA         /**< [out] pi^T A */);

    /** Add columns */
    virtual DSP_RTN_CODE addCols(
    		const double* price,                  /**< [in] price */
    		const double* piA,                    /**< [in] pi^T A */
    		std::vector<int>& indices,            /**< [in] subproblem indices corresponding to cols*/
			std::vector<int>& statuses,           /**< [in] subproblem solution status */
			std::vector<double>& cxs,             /**< [in] solution times original objective coefficients */
			std::vector<double>& objs,            /**< [in] subproblem objective values */
    		std::vector<CoinPackedVector*>& sols, /**< [in] subproblem solutions */
			int& nadded                           /**< [out] number of columns added */);

    /** Calculate Lagrangian bound */
    virtual DSP_RTN_CODE getLagrangianBound(
    		const double* price,       /**< [in] price */
			std::vector<double>& objs, /**< [in] subproblem objective values */
			double& lb                 /**< [out] lower bound */);

    /** update master */
    virtual DSP_RTN_CODE updateModel(
    		const double* price, /**< [in] price */
			double curLb         /**< [in] current lower bound */);

    /** get warm-start information */
    virtual DSP_RTN_CODE getWarmStartInfo(
    		std::vector<double>& sol, /**< [out] current solution */
    		CoinWarmStartBasis*& ws   /**< [out] warmstart basis */);

    /** set warm-start information */
    virtual DSP_RTN_CODE setWarmStartInfo(
    		std::vector<double>& sol, /**< [out] current solution */
    		CoinWarmStartBasis*& ws   /**< [out] warmstart basis */);

    /** termination test */
    virtual bool terminationTest(int nnewcols, int itercnt, double relgap);

    /** termination test after column generation*/
    virtual bool terminationTestColgen(std::vector<int>& statuses);

    /** run heuristics */
    virtual DSP_RTN_CODE heuristics();

    bool useCpxBarrier_;

    int phase_; /**< phase 1 or 2? */
    std::vector<CoinPackedVector*> auxcols_;
    std::vector<int> auxcolindices_;

    /**@name original row bounds of the branching constraints,
     * which is equivalent to the original integer column bounds. */
    //@{
    double* rlbd_branch_;
    double* rubd_branch_;
    //@}

    DwWorker* worker_; /**< subproblem solver */

    std::vector<DwCol*> cols_generated_; /**< columns generated */

public:

    int ncols_orig_; /**< number of columns in the original master */

    int nrows_;        /**< number of rows */
    int nrows_orig_;   /**< number of rows in the original master */
    int nrows_branch_; /**< number of rows representing integer columns in the original master */
    int nrows_conv_;   /**< number of rows representing the convexification constraints */
    std::map<int,int> branch_row_to_col_; /**< maps each branching row to column in the original master */

    /**@name original master problem data */
    CoinPackedMatrix* org_mat_; /**< constraint matrix */
    double* org_clbd_;
    double* org_cubd_;
    double* org_obj_;
    char* org_ctype_;
    double* org_rlbd_;
    double* org_rubd_;

    int itercnt_;
};

#endif /* SRC_SOLVER_DANTZIGWOLFE_DWALGO_H_ */
