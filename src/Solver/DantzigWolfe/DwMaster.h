//
// Created by Kibaek Kim on 8/27/16.
//

#ifndef DSP_DWMASTER_H
#define DSP_DWMASTER_H

#include "Solver/DantzigWolfe/DwAlgo.h"
#include "TreeSearch/DspBranch.h"
#include "Solver/DantzigWolfe/DwWorker.h"
#include "Solver/DantzigWolfe/DwWorkerMpi.h"

class DwMaster : public DwAlgo {
public:

    /** constructor with worker */
    DwMaster(DwWorker* worker): DwAlgo(worker) {}

    /** default destructor */
    virtual ~DwMaster() {}

    /** initialize */
    virtual DSP_RTN_CODE init();

    /** solve */
    virtual DSP_RTN_CODE solve();

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

    /** pre-process for heuristics */
    virtual DSP_RTN_CODE preHeuristic(
    		double*& rlbd,    /**< [out] original row lower bounds */
    		double*& rubd,    /**< [out] original row lower bounds */
    		double*& primsol, /**< [out] original primal solution */
    		double& primobj,  /**< [out] original primal objective */
    		double& dualobj,  /**< [out] original dual objective */
    		int& status       /**< [out] original solution status */);

    /** post-process for heuristics */
    virtual DSP_RTN_CODE postHeuristic(
    		double*& rlbd,    /**< [out] original row lower bounds */
    		double*& rubd,    /**< [out] original row lower bounds */
    		double*& primsol, /**< [out] original primal solution */
    		double& primobj,  /**< [out] original primal objective */
    		double& dualobj,  /**< [out] original dual objective */
    		int& status       /**< [out] original solution status */);

    /** trivial heuristic */
    virtual DSP_RTN_CODE heuristicTrivial();

    /** FP-type heuristic */
    virtual DSP_RTN_CODE heuristicFp(int direction);

    /** Dive heuristic */
    virtual DSP_RTN_CODE heuristicDive();

    /** guts of Dive heuristic */
    virtual DSP_RTN_CODE gutsOfDive(
    		std::vector<CoinTriple<int,int,double> > branchList,
    		int depth);

    /** Solve the master with integrality */
    DSP_RTN_CODE solveMip();

};


#endif //DSP_DWMASTER_H
