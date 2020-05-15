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

	/** solve */
	virtual DSP_RTN_CODE solve();

	/** get best dual objective */
	virtual double getBestDualObjective() {return -bestdualobj_;}

	/** get dual objective */
	virtual double getDualObjective() {return -dualobj_;}

	/** A virtual member to return relative approximate gap */
	virtual double getRelApproxGap() {return fabs(primobj_+bestdualobj_) / (1.0e-10 + fabs(bestdualobj_));}

protected:

	/** This creates a master problem. */
	virtual DSP_RTN_CODE createProblem();

	/** Update the proximal bundle center */
	virtual DSP_RTN_CODE updateCenter(double penalty);

	/** solver master */
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

	/** Calculate Lagrangian bound */
	virtual DSP_RTN_CODE getLagrangianBound(
			std::vector<double>& objs /**< [in] subproblem objective values */);

	/** print iteration information */
	virtual void printIterInfo();

	/** create primal master problem */
	virtual DSP_RTN_CODE createPrimalProblem();

	/** create dual master problem */
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

	/** call external solver for solveMaster() */
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
