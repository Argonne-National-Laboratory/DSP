/*
 * DwModelSmip.h
 *
 *  Created on: Jan 20, 2020
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DANTZIGWOLFE_DWMODELSMIP_H_
#define SRC_SOLVER_DANTZIGWOLFE_DWMODELSMIP_H_

#include "Solver/DantzigWolfe/DwModel.h"
#include "Model/TssModel.h"

class DwModelSmip: public DwModel {
public:
	/** default constructor */
	DwModelSmip();

	/** default constructor with solver */
	DwModelSmip(DecSolver* solver);

	/** default destructor */
	virtual ~DwModelSmip();
				
    /** calculate and return reference solution */
    virtual void getRefSol(std::vector<double>& refsol);

protected:

	/** initialize branching */
	virtual DSP_RTN_CODE initBranch();

	/** initialize heuristic */
	virtual DSP_RTN_CODE initHeuristic();

	/** parse primal solution from the master */
	virtual DSP_RTN_CODE parsePrimSolution();

	/** parse primal solution from the last iteration */
	virtual DSP_RTN_CODE parseLastIterSolution();

    /** calculate and return deviation of the reference solution */
    virtual double getMaxDev(std::vector<double> primsol);

private:

    TssModel* tss_;
};

#endif /* SRC_SOLVER_DANTZIGWOLFE_DWMODELSMIP_H_ */
