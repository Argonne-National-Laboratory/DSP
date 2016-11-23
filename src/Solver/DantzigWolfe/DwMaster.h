//
// Created by Kibaek Kim on 8/27/16.
//

#ifndef DSP_DWMASTER_H
#define DSP_DWMASTER_H

#include "Solver/DantzigWolfe/DwAlgo.h"
#include "Solver/DantzigWolfe/DwFeasPump.h"
#include "TreeSearch/DspBranch.h"

class DwMaster : public DwAlgo {
public:

    /** default constructor */
    DwMaster(
    		DecModel *   model,  /**< model pointer */
            DspParams *  par,    /**< parameters */
            DspMessage * message /**< message pointer */):
	DwAlgo(model, par, message), fpump_(NULL) {}

    /** default destructor */
    virtual ~DwMaster() {
    	FREE_PTR(fpump_);
    }

    /** initialize */
    virtual DSP_RTN_CODE init();

	/** The function chooses branching objects and returns the pointers. */
	virtual bool chooseBranchingObjects(
			DspBranch*& branchingUp, /**< [out] branching-up object */
			DspBranch*& branchingDn  /**< [out] branching-down object */);

	/** set branching objects */
	void setBranchingObjects(const DspBranch* branchobj);

protected:

    /** This creates a master problem. */
	virtual DSP_RTN_CODE createProblem();

    /** Run heuristics */
    virtual DSP_RTN_CODE heuristics();

    /** Solve the master with integrality */
    DSP_RTN_CODE solveMip();

public:

    /** redefine getNumCols() */
    int getNumCols() {return ncols_orig_;}

protected:

    DwFeasPump* fpump_;
};


#endif //DSP_DWMASTER_H
