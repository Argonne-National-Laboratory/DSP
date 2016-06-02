/*
 * DdWorkerLB.h
 *
 *  Created on: Apr 4, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DUALDECOMP_DDWORKERLB_H_
#define SRC_SOLVER_DUALDECOMP_DDWORKERLB_H_

#include <Solver/DualDecomp/DdWorker.h>

class DdWorkerLB: public DdWorker {

	friend class DdMWSync;
	friend class DdMWAsync;

public:

	/** constructor */
	DdWorkerLB(DspParams * par, DecModel * model, DspMessage * message);

	/** destructor */
	virtual ~DdWorkerLB();

	/** initialize */
	virtual DSP_RTN_CODE init();

	/** solve */
	virtual DSP_RTN_CODE solve();

public:

	/** get worker type */
	virtual int getType() {return LB;}

protected:

	/** create problem */
	virtual DSP_RTN_CODE createProblem();

protected:

	vector<DdSub*> subprobs_; /**< set of subproblems */
};

#endif /* SRC_SOLVER_DUALDECOMP_DDWORKERLB_H_ */
