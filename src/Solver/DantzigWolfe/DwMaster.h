//
// Created by Kibaek Kim on 8/27/16.
//

#ifndef DSP_DWMASTER_H
#define DSP_DWMASTER_H

/** Coin */
#include "CoinWarmStartBasis.hpp"
/** Dsp */
#include "TreeSearch/DspBranch.h"
#include "Solver/DecSolver.h"
#include "Solver/DantzigWolfe/DwCol.h"
#include "Solver/DantzigWolfe/DwWorker.h"
#include "Solver/DantzigWolfe/DwWorkerMpi.h"

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
class DwMaster: public DecSolver {
public:

    /** constructor with worker */
    DwMaster(DwWorker* worker);

    /** default destructor */
    virtual ~DwMaster();

    /** initialize */
    virtual DSP_RTN_CODE init();

    /** solve */
    virtual DSP_RTN_CODE solve();

    /** finalize */
    virtual DSP_RTN_CODE finalize() {return DSP_RTN_OK;}

    /** switch to phase 2 */
    DSP_RTN_CODE switchToPhase2();

	/** set branching objects */
	virtual void setBranchingObjects(const DspBranch* branchobj);

	/** set best primal solution */
	void setBestPrimalSolution(const double* solution);

	/** set primal solution */
	void setPrimalSolution(const double* solution);

protected:

	/** Generate initial columns */
	virtual DSP_RTN_CODE initialColumns();

    /** This creates a master problem. */
	virtual DSP_RTN_CODE createProblem();

	/** initial solver setting */
	virtual DSP_RTN_CODE initialOsiSolver();

    /** solve phase 1 */
    virtual DSP_RTN_CODE solvePhase1();

    /** solve phase 2 */
    virtual DSP_RTN_CODE solvePhase2();

    /** guts of solve */
    virtual DSP_RTN_CODE gutsOfSolve();

    /** solver master */
    virtual DSP_RTN_CODE solveMaster();

    /** restore columns: adding all the columns back */
    virtual DSP_RTN_CODE restoreCols();

    /** reduce columns (e.g., reduced cost fixing) */
    virtual DSP_RTN_CODE reduceCols();

    /** generate columns */
    virtual DSP_RTN_CODE generateCols(
			int& ncols /**< [out] number of columns generated */);

    /** calculate piA */
    virtual DSP_RTN_CODE calculatePiA(
			std::vector<double>& piA /**< [out] pi^T A */);

    /** Add columns */
    virtual DSP_RTN_CODE addCols(
    		const double* piA,                    /**< [in] pi^T A */
    		std::vector<int>& indices,            /**< [in] subproblem indices corresponding to cols*/
			std::vector<int>& statuses,           /**< [in] subproblem solution status */
			std::vector<double>& cxs,             /**< [in] solution times original objective coefficients */
			std::vector<double>& objs,            /**< [in] subproblem objective values */
    		std::vector<CoinPackedVector*>& sols, /**< [in] subproblem solutions */
			int& nadded                           /**< [out] number of columns added */);

    /** Calculate Lagrangian bound */
    virtual DSP_RTN_CODE getLagrangianBound(
			std::vector<double>& objs /**< [in] subproblem objective values */);

    /** update master */
    virtual DSP_RTN_CODE updateModel();

    /** termination test */
    virtual bool terminationTest(int nnewcols);
#if 0
    /** Run heuristics */
    virtual DSP_RTN_CODE heuristics();

    /** pre-process for heuristics */
    virtual DSP_RTN_CODE preHeuristic(
    		double*& rlbd,    /**< [out] original row lower bounds */
    		double*& rubd,    /**< [out] original row lower bounds */
    		double*& primsol, /**< [out] original primal solution */
    		double& primobj,  /**< [out] original primal objective */
    		double& dualobj,  /**< [out] original dual objective */
    		int& status       /**< [out] original solution status */);

    /** post-process for heuristics */
    virtual DSP_RTN_CODE postHeuristic(
    		double*& rlbd,    /**< [out] original row lower bounds */
    		double*& rubd,    /**< [out] original row lower bounds */
    		double*& primsol, /**< [out] original primal solution */
    		double& primobj,  /**< [out] original primal objective */
    		double& dualobj,  /**< [out] original dual objective */
    		int& status       /**< [out] original solution status */);

    /** trivial heuristic */
    virtual DSP_RTN_CODE heuristicTrivial();

    /** FP-type heuristic */
    virtual DSP_RTN_CODE heuristicFp(int direction);

    /** Dive heuristic */
    virtual DSP_RTN_CODE heuristicDive();

    /** guts of Dive heuristic */
    virtual DSP_RTN_CODE gutsOfDive(
    		std::vector<CoinTriple<int,int,double> > branchList,
    		int depth);
#endif
    /** print iteration information */
    virtual void printIterInfo();

protected:

    bool useBarrier_;

    int phase_; /**< phase 1 or 2? */
    std::vector<int> auxcolindices_;

    DwWorker* worker_; /**< subproblem solver */

    const double feastol_ = 1.0e-5;

public:

    int ncols_orig_; /**< number of columns in the original master */
    int ncols_start_; /**< index that generated columns start added. */

    int nrows_;        /**< number of rows */
    int nrows_orig_;   /**< number of rows in the original master */
    int nrows_conv_;   /**< number of rows representing the convexification constraints */
    int nrows_core_;   /**< nrows_orig_ + nrows_conv_ */
    int nrows_branch_; /**< number of rows representing integer columns in the original master */
    std::map<int,int> branch_row_to_col_; /**< maps each branching row to column in the original master */

    std::vector<DwCol*> cols_generated_; /**< columns generated */

    /**@name original master problem data */
    CoinPackedMatrix* mat_orig_; /**< constraint matrix */
    double* clbd_orig_;
    double* cubd_orig_;
    double* obj_orig_;
    char* ctype_orig_;
    double* rlbd_orig_;
    double* rubd_orig_;

	double* clbd_node_; /** current column lower bounds */
	double* cubd_node_; /** current column upper bounds */

    int itercnt_;
    int ngenerated_;

    /**@name Time stamps */
    double t_total_; /**< total time */
    double t_master_; /**< master solution time */
    double t_colgen_; /**< column generation time */

protected:

    std::vector<int> status_subs_; /**< subproblem status */

};


#endif //DSP_DWMASTER_H
