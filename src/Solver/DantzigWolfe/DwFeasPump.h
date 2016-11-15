/*
 * DwFeasPump.h
 *
 *  Created on: Oct 27, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DANTZIGWOLFE_DWFEASPUMP_H_
#define SRC_SOLVER_DANTZIGWOLFE_DWFEASPUMP_H_

#include "Solver/DantzigWolfe/DwAlgo.h"

class DwFeasPump: public DwAlgo {
public:
    /** default constructor */
	DwFeasPump(
			DecModel *   model,  /**< model pointer */
            DspParams *  par,    /**< parameters */
            DspMessage * message /**< message pointer */);

    /** default destructor */
	virtual ~DwFeasPump() {
		for (unsigned i = 0; i < visited_.size(); ++i)
			FREE_PTR(visited_[i]);
	}

    /** solve */
    virtual DSP_RTN_CODE solve();

    /** copy columns */
    virtual DSP_RTN_CODE copyColumns(std::vector<DwCol*>& cols);

    /** set branch row bounds */
    virtual DSP_RTN_CODE setBranchRowBounds(const double* rlbd, const double* rubd);

    /** get new columns generated during this heuristic */
    virtual DSP_RTN_CODE getNewCols(std::vector<DwCol*>& cols);

protected:

    /** This creates a master problem. */
    virtual DSP_RTN_CODE createProblem();

    /** set of solutions visited */
	vector<CoinPackedVector*> visited_;
};

#endif /* SRC_SOLVER_DANTZIGWOLFE_DWFEASPUMP_H_ */
