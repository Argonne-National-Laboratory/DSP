/*
 * DwMasterTrLight.h
 *
 *  Created on: Dec 13, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DANTZIGWOLFE_DWMASTERTRLIGHT_H_
#define SRC_SOLVER_DANTZIGWOLFE_DWMASTERTRLIGHT_H_

#include <DantzigWolfe/DwMasterTr.h>

class DwMasterTrLight: public DwMasterTr {
public:
    /** constructor with worker*/
	DwMasterTrLight(DwWorker* worker);

    /** default destructor */
	virtual ~DwMasterTrLight();

	/** The function chooses branching objects and returns the pointers. */
	virtual bool chooseBranchingObjects(
			DspBranch*& branchingUp, /**< [out] branching-up object */
			DspBranch*& branchingDn  /**< [out] branching-down object */);

	/** set branching objects */
	void setBranchingObjects(const DspBranch* branchobj);

protected:

    /** This creates a master problem. */
	virtual DSP_RTN_CODE createProblem();

private:

	double* node_clbd_; /** current column lower bounds */
	double* node_cubd_; /** current column upper bounds */
};

#endif /* SRC_SOLVER_DANTZIGWOLFE_DWMASTERTRLIGHT_H_ */
