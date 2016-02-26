/*
 * DdWorker.h
 *
 *  Created on: Feb 9, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DUALDECOMP_DDWORKER_H_
#define SRC_SOLVER_DUALDECOMP_DDWORKER_H_

/** DSP */
#include "Solver/DualDecomp/DdSolver.h"
#include "Solver/DualDecomp/DdSub.h"

class DdWorker : public DdSolver {
public:

	/** constructor */
	DdWorker(DspParams * par, DecModel * model, StoMessage * message);

	/** destructor */
	virtual ~DdWorker();

	/** initialize */
	virtual STO_RTN_CODE init();

	/** solve */
	virtual STO_RTN_CODE solve();

	/**
	 * Messages to communicate with master
	 *
	 * Structure of the message to receive:
	 *   1 signal to either continue or stop
	 *   [for each subproblem]
	 *   2 theta
	 *   3 lambda
	 *
	 * Structure of the message to send:
	 *   1 number of subproblems
	 *   [for each subproblem]
	 *   2 subproblem index
	 *   3 primal objective
	 *   4 dual objective
	 *   5 coupling column part of the solution
	 */

	/** receive message from worker */
	virtual STO_RTN_CODE recvMessage(int size, double * message);

protected:

	/** create problem */
	virtual STO_RTN_CODE createProblem();

	/** create message */
	virtual STO_RTN_CODE createMessage();

protected:

	vector<DdSub*> subprobs_; /**< set of subproblems */
};

#endif /* SRC_SOLVER_DUALDECOMP_DDWORKER_H_ */
