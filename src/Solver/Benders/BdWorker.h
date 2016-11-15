/*
 * BdWorker.h
 *
 *  Created on: Feb 18, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_BENDERS_BDWORKER_H_
#define SRC_SOLVER_BENDERS_BDWORKER_H_

#include "Solver/Benders/BdSub.h"

class BdWorker {
public:

	/** constructor */
	BdWorker(DecModel * model, DspParams * par, DspMessage * message);

	/** destructor */
	virtual ~BdWorker();

	/** generate cuts */
	DSP_RTN_CODE generateCuts(int nx, int naux, const double* x, OsiCuts& cs);

protected:

	/** aggregate cuts, if necessary, based on the number of auxiliary variables
	 * introduced to the master.
	 *
	 * This may be derived to communicate cuts in MPI parallel implementation.
	 */
	virtual DSP_RTN_CODE collectCuts(int nx, int naux, double** cut, double* rhs, OsiCuts& cs);

protected:

	DecModel * model_;     /**< DecModel object */
	DspParams * par_;      /**< parameters */
	DspMessage * message_; /**< message */

	BdSub * bdsub_; /**< benders subproblem */

	int parProcIdxSize_; /**< number of subproblems for this worker */
	int * parProcIdx_;   /**< subproblem indices for this worker */
};

#endif /* SRC_SOLVER_BENDERS_BDWORKER_H_ */
