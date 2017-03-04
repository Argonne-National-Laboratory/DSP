/*
 * DwModel.h
 *
 *  Created on: Feb 13, 2017
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DANTZIGWOLFE_DWMODEL_H_
#define SRC_SOLVER_DANTZIGWOLFE_DWMODEL_H_

#include <TreeSearch/DspModel.h>
#include <DantzigWolfe/DwMaster.h>

class DwModel: public DspModel {
public:
	/** default constructor */
	DwModel();

	/** default constructor with solver */
	DwModel(DecSolver* solver);

	/** default destructor */
	virtual ~DwModel();

	/** solve model */
    virtual DSP_RTN_CODE solve();

    virtual bool chooseBranchingObjects(
    			DspBranch*& branchingUp, /**< [out] branching-up object */
    			DspBranch*& branchingDn  /**< [out] branching-down object */);

private:

    DwMaster* master_;
    double infeasibility_;
};

#endif /* SRC_SOLVER_DANTZIGWOLFE_DWMODEL_H_ */
