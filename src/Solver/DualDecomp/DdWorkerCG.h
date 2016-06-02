/*
 * DdWorkerCG.h
 *
 *  Created on: May 23, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DUALDECOMP_DDWORKERCG_H_
#define SRC_SOLVER_DUALDECOMP_DDWORKERCG_H_

#include "Solver/DualDecomp/DdWorker.h"

class DdWorkerCG: public DdWorker {
public:
	/** constructor */
	DdWorkerCG(DspParams * par, DecModel * model, DspMessage * message);

	/** destructor */
	virtual ~DdWorkerCG();

protected:
	int parProcIdxSize_; /**< number of subproblems for this worker */
	int * parProcIdx_;   /**< subproblem indices for this worker */
};

#endif /* SRC_SOLVER_DUALDECOMP_DDWORKERCG_H_ */
