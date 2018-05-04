/*
 * DwBranchInt.h
 *
 *  Created on: May 1, 2018
 *      Author: Kibaek Kim
 */

#ifndef SRC_SOLVER_DANTZIGWOLFE_DWBRANCHINT_H_
#define SRC_SOLVER_DANTZIGWOLFE_DWBRANCHINT_H_

#include <DantzigWolfe/DwBranch.h>

class DwBranchInt : public DwBranch {
public:
	/** default constructor */
	DwBranchInt() : DwBranch() {}

	/** default constructor with solver */
	DwBranchInt(DwModel* model) : DwBranch(model) {}

	/** default destructor */
	virtual ~DwBranchInt() {
	}

    virtual bool chooseBranchingObjects(
    			std::vector<DspBranchObj*>& branchingObjs /**< [out] branching objects */);

};

#endif /* SRC_SOLVER_DANTZIGWOLFE_DWBRANCH_H_ */
