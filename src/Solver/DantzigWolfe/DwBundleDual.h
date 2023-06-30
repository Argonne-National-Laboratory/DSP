/*
 * DwBundleDual.h
 *
 *  Created on: Feb 20, 2017
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DANTZIGWOLFE_DWBUNDLEDUAL_H_
#define SRC_SOLVER_DANTZIGWOLFE_DWBUNDLEDUAL_H_

#include <memory>
#include <algorithm>
#include <numeric>
#include "Solver/DantzigWolfe/DwMaster.h"

/**
 * @brief This implements a bundle method for solving the master problem.
 * 
 */
class DwBundleDual: public DwMaster {
public:
	/** constructor with worker */
	DwBundleDual(DwWorker* worker);

	/** copy constructor */
	DwBundleDual(const DwBundleDual& rhs);

	/** copy operator */
	DwBundleDual& operator=(const DwBundleDual& rhs);

	virtual DwBundleDual* clone() const {
		return new DwBundleDual(*this);
	}

	/** default destructor */
	virtual ~DwBundleDual();

	/**
	 * @brief This function solves the master problem by using a bundle method.
	 * 
	 * This initializes the best dual solution, the dual solution, and the proximal center in the
	 * quadratic objective function of the dual problem.
	 * 
	 * This calls DwMaster::generateCols() to generate a set of initial columns. These are added as columns
	 * to the primal problem, but added as rows to the dual problem. DwMaster::generateCols() calls 
	 * addCols(), the implementation of which (i.e., DwBundleDual::addCols()) in this class calls 
	 * DwBundleDual::addRows(), where the columns are added to the primal problem object and the rows are
	 * added to the dual proble object.
	 * 
	 * After adding the initial rows/columns, DwMaster::gutsOfSolver() is called.
	 * 
	 * @return DSP_RTN_CODE 
	 */
	virtual DSP_RTN_CODE solve();

	/** get best dual objective */
	virtual double getBestDualObjective() {return -bestdualobj_;}

	/** get dual objective */
	virtual double getDualObjective() {return -dualobj_;}

	/** A virtual member to return relative approximate gap */
	virtual double getRelApproxGap() {return fabs(primobj_+bestdualobj_) / (1.0e-10 + fabs(bestdualobj_));}

protected:

	/**
	 * @brief This overwrites DwMaster::createProblem() and creates the primal and dual problem objects. 
	 * 
	 * This creates the primal and dual problem objects and sets the column bounds for the current node
	 * subproblem. phase_ is always set to 2, because the Bundle method does not consider two phases.
	 * The best dual objective value is initialized to COIN_DBL_MAX;
	 * 
	 * @return DSP_RTN_CODE 
	 */
	virtual DSP_RTN_CODE createProblem();

	/** Update the proximal bundle center */
	virtual DSP_RTN_CODE updateCenter(double penalty);

	/**
	 * @brief This solves the master problems in both primal and dual.
	 * 
	 * This function is called from DwMaster::gutsOfSolve(). The dual problem is solved first. If the solution 
	 * status is either optimal or feasible, this function updates the bundle parameters and solves the primal 
	 * problem.
	 * 
	 * @return DSP_RTN_CODE 
	 */
	virtual DSP_RTN_CODE solveMaster();

	/** update master */
	virtual DSP_RTN_CODE updateModel();

	/** termination test */
	virtual bool terminationTest();

    /** reduce columns */
    virtual DSP_RTN_CODE reduceCols(int &num_removed);

    /** restore columns: adding all the columns back */
    virtual DSP_RTN_CODE restoreCols(int &num_restored);

	/** reduce columns (e.g., reduced cost fixing) */
	virtual DSP_RTN_CODE reduceCols() {return DSP_RTN_OK;}

	/** Add columns */
	virtual DSP_RTN_CODE addCols(
			std::vector<int>& indices,           /**< [in] subproblem indices corresponding to cols*/
			std::vector<int>& statuses,          /**< [in] subproblem solution status */
			std::vector<double>& cxs,            /**< [in] solution times original objective coefficients */
			std::vector<double>& objs,           /**< [in] subproblem objective values */
			std::vector<CoinPackedVector*>& sols /**< [in] subproblem solutions */);

	/** Add rows */
	virtual DSP_RTN_CODE addRows(
			std::vector<int>& indices,           /**< [in] subproblem indices corresponding to cols*/
			std::vector<int>& statuses,          /**< [in] subproblem solution status */
			std::vector<double>& cxs,            /**< [in] solution times original objective coefficients */
			std::vector<double>& objs,           /**< [in] subproblem objective values */
			std::vector<CoinPackedVector*>& sols /**< [in] subproblem solutions */);

	/**
	 * @brief This computes the Lagrangian dual bound.
	 * 
	 * The sign of the Lagrangian dual bound is same as the original objective function.
	 */
	virtual DSP_RTN_CODE getLagrangianBound(
			std::vector<double>& objs /**< [in] subproblem objective values */);

	/** print iteration information */
	virtual void printIterInfo();

	/**
	 * @brief This creates a primal (minimization) problem object.
	 * 
	 * The primal problem is the restricted master problem of Dantzig-Wolfe decomposition.
	 * The problem consists of two sets of constraints.
	 * The first set of constraints are reserved (empty initially) to represents the convex 
	 * combination of columns generated in this method. The other constraints are reserved to 
	 * represent the original constraints of the master block. The solver interface is created.
	 * 
	 * @return DSP_RTN_CODE 
	 */
	virtual DSP_RTN_CODE createPrimalProblem();

	/**
	 * @brief This creates a dual (maximization) problem object.
	 * 
	 * The problem is initialized with (i) the auxiliary variables to represent the upper approximation
	 * of the blocks and (ii) the variables to represent the Lagrangian multipliers of the master block
	 * constraints. The quadratic objective function is created for the proximal term with weight u_;
	 * 
	 * The problem is minimized with the negated objective function.
	 * 
	 * @return DSP_RTN_CODE 
	 */
	virtual DSP_RTN_CODE createDualProblem();

	typedef std::pair<int,double> pairIntDbl;
	static bool compPair ( const pairIntDbl& l, const pairIntDbl& r) {
		return l.second < r.second;
	}

	const double mL_ = 0.1; /**< parameter for serious step, = (0,0.5) */
	const double mR_ = 0.5; /**< parameter for deciding whether u_ is too large, = (mL_,1)*/
	const double umin_ = 1.0e-4; /**< minimum weight */
	const double umax_ = 1.0e+4; /**< maximum weight */
	double v_; /**< predicted ascent */
	int counter_; /**< number of serious steps or null steps since the latest change of u_ */
	double u_; /**< bundle weight */
	double eps_; /** variation estimate */
	std::vector<double> d_; /**< direction = dual - center */
	std::vector<double> p_; /**< aggregate subgradient */
	double absp_; /**< norm of p_ */
	double alpha_; /**< aggregate linearization error at the best dual */
	double linerr_; /**< linearization error */

	double prev_dualobj_; /**< dual objective at the previous iteration */
	int nstalls_; /**< number of iterations making no progress on objective value */

	int numFixedRows_; /**< number of fixed rows in the dual master */

	std::shared_ptr<DspOsi> primal_si_;

	//@{
	/** functions specific to external solver */

	/** initialize dual solver */
	virtual void initDualSolver(
			const CoinPackedMatrix& m, 
			std::vector<double>& clbd, 
			std::vector<double>& cubd, 
			std::vector<double>& obj, 
			std::vector<double>& rlbd, 
			std::vector<double>& rubd);

	/**
	 * @brief This solves the dual problem and returns a status.
	 * 
	 * If CPLEX is used as external solver, extra effort is made when the problem suffers from
	 * issues in numerical stability.
	 * 
	 * @return DSP_RTN_CODE 
	 */
	virtual DSP_RTN_CODE callMasterSolver();

	/** assign master solution */
	virtual void assignMasterSolution(std::vector<double>& sol);

	/** get objective value */
	virtual double getObjValue() { return osi_->si_->getObjValue(); }

	/** add row to the dual master */
	virtual void addDualRow(const CoinPackedVector& v, const double lb, const double ub) {
		osi_->si_->addRow(v, lb, ub);
	}

    /** remove all columns in the DW master */
    virtual void removeAllCols();

	/** remove all columns in the primal master */
	virtual void removeAllPrimCols();

	/** remove all rows in the dual master */
	virtual void removeAllDualRows();

	/** remove branching rows from the DW master */
	virtual void removeBranchingRows();

	/** add branching rows and columns */
	//virtual void addBranchingRowsCols(const DspBranchObj* branchobj);

    /** add branching row */
    virtual void addBranchingRow(double lb, double ub);

    /** add branching column */
    virtual void addBranchingCol(const CoinPackedVector& col, double obj);

	//@}
};

#endif /* SRC_SOLVER_DANTZIGWOLFE_DWBUNDLEDUAL_H_ */
