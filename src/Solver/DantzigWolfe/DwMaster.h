//
// Created by Kibaek Kim on 8/27/16.
//

#ifndef DSP_DWMASTER_H
#define DSP_DWMASTER_H

//#define USE_ROW_TO_COL

#include "CoinWarmStartBasis.hpp"
#include "Solver/DecSolver.h"
#include "Solver/DantzigWolfe/DwCol.h"
#include "Solver/DantzigWolfe/DwWorker.h"
#ifdef DSP_HAS_MPI
#include "Solver/DantzigWolfe/DwWorkerMpi.h"
#endif /* DSP_HAS_MPI */

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

    /** copy constructor */
    DwMaster(const DwMaster& rhs);

    /** copy operator */
    DwMaster& operator=(const DwMaster& rhs);

	virtual DwMaster* clone() const {
		return new DwMaster(*this);
	}

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

    /** get best primal solution (in original space) found during iteration */
    virtual const double * getBestPrimalSolutionOrig() {return &bestprimsol_orig_[0];}

    /** get best primal solution (in original space) found during iteration */
    virtual std::vector<CoinPackedVector*>& getLastSubprobSolutions() {return recent_subsols_;}

	/** set branching objects */
	virtual void setBranchingObjects(const DspBranchObj* branchobj);

protected:

    /** create DspOsi for master */
    virtual DspOsi * createDspOsi();

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

    /** update column vector */
    virtual DSP_RTN_CODE updateCol(DwCol* col);

    /** restore columns: adding all the columns back */
    virtual DSP_RTN_CODE restoreCols(int &num_restored);

    /** reduce columns (e.g., reduced cost fixing) */
    virtual DSP_RTN_CODE reduceCols(int &num_removed);

    /** generate columns */
    virtual DSP_RTN_CODE generateCols();

    /** generate columns by fixing first-stage variables of SMIP */
    virtual DSP_RTN_CODE generateColsByFix(
        int nsols /**< [in] number of solutions to evaluate in LIFO way */);

    /** calculate piA */
    virtual DSP_RTN_CODE calculatePiA(
			std::vector<double>& piA /**< [out] pi^T A */);

    /** Add columns */
    virtual DSP_RTN_CODE addCols(
    		std::vector<int>& indices,           /**< [in] subproblem indices corresponding to cols*/
			std::vector<int>& statuses,          /**< [in] subproblem solution status */
			std::vector<double>& cxs,            /**< [in] solution times original objective coefficients */
			std::vector<double>& objs,           /**< [in] subproblem objective values */
    		std::vector<CoinPackedVector*>& sols /**< [in] subproblem solutions */);

    /** remove all columns in the master */
    virtual void removeAllCols();

    /** remove branching rows */
    virtual void removeBranchingRows();

    /** add branching rows and columns */
    virtual void addBranchingRowsCols(const DspBranchObj* branchobj);

    /** add branching row */
    virtual void addBranchingRow(double lb, double ub);

    /** add branching column */
    virtual void addBranchingCol(const CoinPackedVector& col, double obj);

    /** Calculate Lagrangian bound */
    virtual DSP_RTN_CODE getLagrangianBound(
			std::vector<double>& objs /**< [in] subproblem objective values */);

    /** update master */
    virtual DSP_RTN_CODE updateModel();

    /** termination test */
    virtual bool terminationTest();
    
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
#ifdef USE_ROW_TO_COL
    std::map<int,int> branch_row_to_col_; /**< maps each branching row to column in the original master */
#else
    std::map<int,CoinPackedVector> branch_row_to_vec_;
#endif

    std::vector<DwCol*> cols_generated_; /**< columns generated */

    /**@name original master problem data */
    CoinPackedMatrix* mat_orig_; /**< constraint matrix */
    std::vector<double> clbd_orig_;
    std::vector<double> cubd_orig_;
    std::vector<double> obj_orig_;
    std::vector<char> ctype_orig_;
    std::vector<double> rlbd_orig_;
    std::vector<double> rubd_orig_;


	std::vector<double> clbd_node_; /** current column lower bounds */
    std::vector<double> cubd_node_; /** current column upper bounds */

    const DspBranchObj* branchObj_; /** branching object (not owned by this class) */


    std::vector<double> bestprimsol_orig_; /** best primal solution in original space */
    std::vector<CoinPackedVector*> recent_subsols_;

    int itercnt_;
    int ngenerated_;

    std::vector<double> log_time_;
    std::vector<double> log_bestdual_bounds_;
    std::vector<double> log_bestprim_bounds_;

    /**@name Time stamps */
    double t_start_; /**< solution start time */
    double t_total_; /**< total time */
    double t_master_; /**< master solution time */
    double t_colgen_; /**< column generation time */

    std::vector<int> status_subs_; /**< subproblem status */
    std::vector<CoinPackedVector*> stored_solutions_; /**< first-stage solutions that are evaluated for upper bounds */

};


#endif //DSP_DWMASTER_H
