/*
 * BdWorker.h
 *
 *  Created on: Feb 18, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_BENDERS_BDWORKER_H_
#define SRC_SOLVER_BENDERS_BDWORKER_H_

#include "Solver/DecSolver.h"
#include "Solver/Benders/BdSub.h"

class BdWorker: public DecSolver {
public:

	/** constructor */
	BdWorker(DspParams * par, DecModel * model, DspMessage * message);

	/** destructor */
	virtual ~BdWorker();

	/** initialize */
	virtual DSP_RTN_CODE init();

	/** solve */
	virtual DSP_RTN_CODE solve();

public:

	/** get BdSub pointer */
	virtual BdSub * getBdSubPtr() {return bdsub_;}

protected:

	/** create problem */
	virtual DSP_RTN_CODE createProblem();

private:

	BdSub * bdsub_; /**< benders subproblem */

	int parProcIdxSize_; /**< number of subproblems for this worker */
	int * parProcIdx_;   /**< subproblem indices for this worker */
};

#endif /* SRC_SOLVER_BENDERS_BDWORKER_H_ */
