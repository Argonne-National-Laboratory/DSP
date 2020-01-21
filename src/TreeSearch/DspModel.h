/*
 * DspModel.h
 *
 *  Created on: Aug 31, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_TREESEARCH_DSPMODEL_H_
#define SRC_TREESEARCH_DSPMODEL_H_

/** Coin */
#include "AlpsModel.h"
/** Dsp */
#include "Solver/DecSolver.h"
#include "TreeSearch/DspHeuristic.h"

/**
 * This implements a model class for Coin-Alps library. This class provides a wrapper
 * for model class and data specific to a decomposition method.
 */
class DspModel: public AlpsModel {
public:

	/** default constructor */
	DspModel();

	/** constructor with solver */
	DspModel(DecSolver* solver);

	/** default destructor */
	virtual ~DspModel();

	/** initialize */
    virtual DSP_RTN_CODE init() {return DSP_RTN_OK;}

	/** solve model */
    virtual DSP_RTN_CODE solve();

    /** Create the root node. Default: do nothing */
    virtual AlpsTreeNode * createRoot();

    /** Return true if all nodes on this process can be fathomed.*/
    virtual bool fathomAllNodes();

public:

    /** returns a pointer to the solver object */
    DecSolver* getSolver() {return solver_;}

    DspParams* getParPtr() {return par_;}
    int getStatus() {return status_;}
    double getBestPrimalObjective() {return bestprimobj_;}
    double getPrimalObjective() {return primobj_;}
    double getBestDualObjective() {return bestdualobj_;}
    double getDualObjective() {return dualobj_;}
    std::vector<double>& getBestPrimalSolution() {return bestprimsol_;}
    std::vector<double>& getPrimalSolution() {return primsol_;}
    double infeasibility() {return infeasibility_;}
    double feastol() {return feastol_;}

    virtual bool chooseBranchingObjects(
    			std::vector<DspBranchObj*>& branchingObjs /**< [out] branching objects */) {
    	return false;
    }

    void setIterLimit(int n) {solver_->setIterLimit(n);}
    void setTimeLimit(double t) {solver_->setTimeLimit(t);}
    void setBestPrimalObjective(double val) {bestprimobj_=val;}
    void setBestDualObjective(double val) {bestdualobj_=val;}
    void setBranchingObjects(const DspBranchObj* branchobj) {
    	solver_->setBranchingObjects(branchobj);
    }

protected:

    DecSolver* solver_; /**< decomposition solver */
    DspParams* par_;

    std::vector<DspHeuristic*> heuristics_;
#if 0
    int ncols_orig_;
    std::vector<char> ctype_orig_;
    std::vector<double> clbd_orig_;
    std::vector<double> cubd_orig_;

    std::vector<double> clbd_node_;
    std::vector<double> cubd_node_;

    std::vector<double> primsol_;
#endif
    int status_;
    double bestprimobj_;
    double bestdualobj_;
    double primobj_;
    double dualobj_;
    std::vector<double> bestprimsol_; /**< integer feasible primal solution */
    std::vector<double> primsol_; /**< primal solution (may not be integer feasible) */
    double infeasibility_;
    double feastol_ = 1.e-6; /**< feasibility tolerance */
};

#endif /* SRC_TREESEARCH_DSPMODEL_H_ */
