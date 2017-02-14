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
    int getStatus() {return solver_->getStatus();}
    int getNumCols() {return solver_->getNumCols();}
    double getBestPrimalObjective() {return solver_->getBestPrimalObjective();}
    double getPrimalObjective() {return solver_->getPrimalObjective();}
    double getDualObjective() {return solver_->getDualObjective();}
    const double* getBestPrimalSolution() {return solver_->getBestPrimalSolution();}
    const double* getPrimalSolution() {return solver_->getPrimalSolution();}

    bool chooseBranchingObjects(
    			DspBranch*& branchingUp, /**< [out] branching-up object */
    			DspBranch*& branchingDn  /**< [out] branching-down object */);

    void setIterLimit(int n) {solver_->setIterLimit(n);}
    void setTimeLimit(double t) {solver_->setTimeLimit(t);}
    void setBranchingObjects(const DspBranch* branchobj) {
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

    double bestprimobj_;
    double bestdualobj_;
    std::vector<double> primsol_;
#endif
};

#endif /* SRC_TREESEARCH_DSPMODEL_H_ */
