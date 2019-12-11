/*
 * DwBranchGenDisj.h
 *
 *  Created on: June 14, 2018
 *      Author: Kibaek Kim
 * 
 *  This implements the branching on general disjunctions.
 */

#ifndef SRC_SOLVER_DANTZIGWOLFE_DWBRANCHGENDISJ_H_
#define SRC_SOLVER_DANTZIGWOLFE_DWBRANCHGENDISJ_H_

#include "Solver/DantzigWolfe/DwBranch.h"
#include "Model/TssModel.h"

class DwBranchGenDisj : public DwBranch {
public:
	/** default constructor */
	DwBranchGenDisj() : DwBranch(), tss_(NULL), master_(NULL) {}

	/** default constructor with solver */
	DwBranchGenDisj(DwModel* model);

	/** default destructor */
	virtual ~DwBranchGenDisj() {
		tss_ = NULL;
		master_ = NULL;
	}

    virtual bool chooseBranchingObjects(
    			std::vector<DspBranchObj*>& branchingObjs /**< [out] branching objects */);

    /** calculate and return separating hyperplane */
    virtual void getSepHyperplane(CoinPackedVector** vec, double& rhs);

protected:

	/** epsilon value for branching on continuous variables */
	double epsilon_ = 1.0e-6;

	TssModel* tss_;
	DwMaster* master_;
};

#endif /* SRC_SOLVER_DANTZIGWOLFE_DWBRANCHGENDISJ_H_ */
