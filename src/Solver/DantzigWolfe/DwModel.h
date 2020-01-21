/*
 * DwModel.h
 *
 *  Created on: Feb 13, 2017
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DANTZIGWOLFE_DWMODEL_H_
#define SRC_SOLVER_DANTZIGWOLFE_DWMODEL_H_

#include "TreeSearch/DspModel.h"
#include "Solver/DantzigWolfe/DwMaster.h"

class DwBranch;

class DwModel: public DspModel {
public:
	/** default constructor */
	DwModel();

	/** default constructor with solver */
	DwModel(DecSolver* solver);

	/** default destructor */
	virtual ~DwModel();

	/** initialize */
	virtual DSP_RTN_CODE init();

	/** solve model */
    virtual DSP_RTN_CODE solve();

    virtual bool chooseBranchingObjects(
    			std::vector<DspBranchObj*>& branchingObjs /**< [out] branching objects */);
				
    /** calculate and return reference solution */
    virtual void getRefSol(std::vector<double>& refsol);

protected:

	/** initialize branching */
	virtual DSP_RTN_CODE initBranch();

	/** initialize heuristic */
	virtual DSP_RTN_CODE initHeuristic();

	/** parse primal solution from the master */
	virtual DSP_RTN_CODE parsePrimSolution();

	/** parse Dantzig-Wolfe solution from the master */
	virtual DSP_RTN_CODE parseDantzigWolfeSolution();

public:

	double heuristic_time_elapsed_;

protected:

    DwBranch* branch_;
};

#endif /* SRC_SOLVER_DANTZIGWOLFE_DWMODEL_H_ */
