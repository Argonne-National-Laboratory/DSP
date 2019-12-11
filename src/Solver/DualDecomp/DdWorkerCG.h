/*
 * DdWorkerCG.h
 *
 *  Created on: May 23, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DUALDECOMP_DDWORKERCG_H_
#define SRC_SOLVER_DUALDECOMP_DDWORKERCG_H_

#include "Solver/DualDecomp/DdWorker.h"

/** An abstract cut generation class */
class DdWorkerCG: public DdWorker {
public:
	/** A default constructor. */
	DdWorkerCG(
			DecModel *   model,  /**< model pointer */
			DspParams *  par,    /**< parameter pointer */
			DspMessage * message /**< message pointer */);

	/** A copy constrcutor. */
	DdWorkerCG(const DdWorkerCG& rhs);

	/** A default destructor. */
	virtual ~DdWorkerCG() {}

	/** A clone function */
	virtual DdWorkerCG* clone() const {
		return new DdWorkerCG(*this);
	}

	/** A virtual member to return the worker type. */
	virtual int getType() {return CG;}

protected:
	int parProcIdxSize_; /**< number of subproblems for this worker */
	const int * parProcIdx_;   /**< subproblem indices for this worker */
};

#endif /* SRC_SOLVER_DUALDECOMP_DDWORKERCG_H_ */
