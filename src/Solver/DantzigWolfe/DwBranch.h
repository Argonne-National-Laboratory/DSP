/*
 * DwBranch.h
 *
 *  Created on: May 1, 2018
 *      Author: Kibaek Kim
 */

#ifndef SRC_SOLVER_DANTZIGWOLFE_DWBRANCH_H_
#define SRC_SOLVER_DANTZIGWOLFE_DWBRANCH_H_

#include <DantzigWolfe/DwModel.h>

class DwBranch {
public:
	/** default constructor */
	DwBranch() : model_(NULL) {}

	/** default constructor with solver */
	DwBranch(DwModel* model) : model_(model) {}

	/** default destructor */
	virtual ~DwBranch() {
		model_ = NULL;
	}

	virtual void setModel(DwModel* model) {model_ = model;}

    virtual bool chooseBranchingObjects(
    			std::vector<DspBranchObj*>& branchingObjs /**< [out] branching objects */) = 0;

protected:

    DwModel* model_;
};

#endif /* SRC_SOLVER_DANTZIGWOLFE_DWBRANCH_H_ */
