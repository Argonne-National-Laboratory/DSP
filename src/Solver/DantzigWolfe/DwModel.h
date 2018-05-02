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

class DwBranch;

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
    			std::vector<DspBranchObj*>& branchingObjs /**< [out] branching objects */);

    DwMaster* getMasterPtr() {return master_;}

private:

    DwMaster* master_;
    DwBranch* branch_;
    double infeasibility_;
};

#endif /* SRC_SOLVER_DANTZIGWOLFE_DWMODEL_H_ */
