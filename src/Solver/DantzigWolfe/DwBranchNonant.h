/*
 * DwBranchNonant.h
 *
 *  Created on: May 2, 2018
 *      Author: Kibaek Kim
 * 
 *  This implements the branching on the nonanticipativity constraints, originally proposed in Caroe and Schultz (1999).
 */

#ifndef SRC_SOLVER_DANTZIGWOLFE_DWBRANCHNONANT_H_
#define SRC_SOLVER_DANTZIGWOLFE_DWBRANCHNONANT_H_

#include <DantzigWolfe/DwBranch.h>
#include <Model/TssModel.h>

class DwBranchNonant : public DwBranch {
public:
	/** default constructor */
	DwBranchNonant() : DwBranch(), tss_(NULL), master_(NULL) {}

	/** default constructor with solver */
	DwBranchNonant(DwModel* model);

	/** default destructor */
	virtual ~DwBranchNonant() {
		tss_ = NULL;
		master_ = NULL;
	}

    virtual bool chooseBranchingObjects(
    			std::vector<DspBranchObj*>& branchingObjs /**< [out] branching objects */);

    /** calculate and return reference solution */
    virtual void getRefSol(std::vector<double>& refsol);

    /** calculate and return deviation of the reference solution */
    virtual void getDevSol(std::vector<double>& refsol, std::vector<double>& devsol);

protected:

	/** epsilon value for branching on continuous variables */
	double epsilon_ = 1.0e-6;
	double epsilonBB_ = 1.0e-6;

	TssModel* tss_;
	DwMaster* master_;
};

#endif /* SRC_SOLVER_DANTZIGWOLFE_DWBRANCHNONANT_H_ */
