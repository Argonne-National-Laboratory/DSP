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
#include <DantzigWolfe/DwMaster.h>

class DwBundleDual: public DwMaster {
public:
    /** constructor with worker */
	DwBundleDual(DwWorker* worker);

    /** default destructor */
	virtual ~DwBundleDual();

    /** solve */
    virtual DSP_RTN_CODE solve();

	/** set branching objects */
	virtual void setBranchingObjects(const DspBranch* branchobj);

protected:

    /** This creates a master problem. */
	virtual DSP_RTN_CODE createProblem();

	/** Update the proximal bundle center */
	DSP_RTN_CODE updateCenter(double penalty);

    /** solver master */
    virtual DSP_RTN_CODE solveMaster();

    /** update master */
    virtual DSP_RTN_CODE updateModel();

    /** termination test */
    virtual bool terminationTest(int nnewcols);

    /** reduce columns (e.g., reduced cost fixing) */
    virtual DSP_RTN_CODE reduceCols() {return DSP_RTN_OK;}

    /** Add columns */
    virtual DSP_RTN_CODE addCols(
    		const double* piA,                    /**< [in] pi^T A */
    		std::vector<int>& indices,            /**< [in] subproblem indices corresponding to cols*/
			std::vector<int>& statuses,           /**< [in] subproblem solution status */
			std::vector<double>& cxs,             /**< [in] solution times original objective coefficients */
			std::vector<double>& objs,            /**< [in] subproblem objective values */
    		std::vector<CoinPackedVector*>& sols, /**< [in] subproblem solutions */
			int& nadded                           /**< [out] number of columns added */);

    /** Add rows */
    virtual DSP_RTN_CODE addRows(
    		std::vector<int>& indices,            /**< [in] subproblem indices corresponding to cols*/
			std::vector<int>& statuses,           /**< [in] subproblem solution status */
			std::vector<double>& cxs,             /**< [in] solution times original objective coefficients */
			std::vector<double>& objs,            /**< [in] subproblem objective values */
    		std::vector<CoinPackedVector*>& sols, /**< [in] subproblem solutions */
			int& nadded                           /**< [out] number of columns added */);

    /** Calculate Lagrangian bound */
    virtual DSP_RTN_CODE getLagrangianBound(
			std::vector<double>& objs /**< [in] subproblem objective values */);

    /** print iteration information */
    virtual void printIterInfo();

private:

    /** create primal master problem */
    DSP_RTN_CODE createPrimalProblem();

    /** create dual master problem */
    DSP_RTN_CODE createDualProblem();

    typedef std::pair<int,double> pairIntDbl;
    static bool compPair ( const pairIntDbl& l, const pairIntDbl& r) {
    	return l.second < r.second;
    }

private:

    const double mL_ = 0.1; /**< parameter for serious step, = (0,0.5) */
    const double mR_ = 0.5; /**< parameter for deciding whether u_ is too large, = (mL_,1)*/
    const double umin_ = 1.0e-10; /**< minimum weight */
    double v_; /**< predicted ascent */
    int counter_; /**< number of serious steps or null steps since the latest change of u_ */
    double u_; /**< bundle weight */
    double eps_; /** variation estimate */
    std::vector<double> p_; /**< aggregate subgradient */
    double absp_; /**< norm of p_ */
    double alpha_; /**< aggregate linearization error at the best dual */
    double linerr_; /**< linearization error */

    std::shared_ptr<OsiSolverInterface> primal_si_;
};

#endif /* SRC_SOLVER_DANTZIGWOLFE_DWBUNDLEDUAL_H_ */
