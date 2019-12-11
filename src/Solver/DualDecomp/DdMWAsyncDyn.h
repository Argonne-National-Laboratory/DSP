/*
 * DdMWAsyncDyn.h
 *
 *  Created on: Feb 20, 2018
 *      Author: Kibaek Kim
 */

#ifndef SRC_SOLVER_DUALDECOMP_DDMWASYNCDYN_H_
#define SRC_SOLVER_DUALDECOMP_DDMWASYNCDYN_H_

#include <Solver/DualDecomp/DdMWAsync.h>

/** A master-worker class for asynchronous dual decomposition with dynamic allocation */
class DdMWAsyncDyn: public DdMWAsync {

public:

	/** A default constructor. */
	DdMWAsyncDyn(
			MPI_Comm     comm,   /**< MPI communicator */
			DecModel *   model,  /**< model pointer */
			DspParams *  par,    /**< parameters */
			DspMessage * message /**< message pointer */);

	/** A copy constructor. */
	DdMWAsyncDyn(const DdMWAsyncDyn& rhs);

	/** A default destructor. */
	virtual ~DdMWAsyncDyn();

	/** A clone function */
	virtual DdMWAsyncDyn* clone() const {
		return new DdMWAsyncDyn(*this);
	}

protected:

	/** choose queue element for evaluating dual variables */
	virtual bool chooseQueueElement(int& qid, double*& qsol, int& nsubprobs, int*& subindex);

	/** set lower bounding workers */
	virtual DSP_RTN_CODE setWorkerLb(DdWorkerLB* workerlb, int nsubprobs, int* subindex, double* buf, double bestprimobj);
};

#endif /* SRC_SOLVER_DUALDECOMP_DDMWASYNCDYN_H_ */
