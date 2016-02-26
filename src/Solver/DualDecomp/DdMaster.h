/*
 * DdMaster.h
 *
 *  Created on: Feb 9, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DUALDECOMP_DDMASTER_H_
#define SRC_SOLVER_DUALDECOMP_DDMASTER_H_

#include "Solver/DualDecomp/DdSolver.h"
#include "SolverInterface/SolverInterface.h"

class DdMaster: public DdSolver {
public:

	/** constructor */
	DdMaster(DspParams * par, DecModel * model, StoMessage * message, int nworkers, int maxnumsubprobs);

	/** destructor */
	virtual ~DdMaster();

	/** initialize */
	virtual STO_RTN_CODE init();

	/**
	 * Messages to communicate with worker
	 *
	 * Structure of the message to send:
	 *   1 signal to either continue or stop
	 *   [for each subproblem]
	 *   2 theta
	 *   3 lambda
	 *
	 * Structure of the message to receive:
	 *   1 number of subproblems
	 *   [for each subproblem]
	 *   2 subproblem index
	 *   3 primal objective
	 *   4 dual objective
	 *   5 coupling column part of the solution
	 */

	/** receive message from worker */
	virtual STO_RTN_CODE recvMessage(int source, int size, double * message);

	/** update problem */
	virtual STO_RTN_CODE updateProblem() = 0;

protected:

	/** create message */
	virtual STO_RTN_CODE createMessage();

public:

	SolverInterface * getSiPtr() {return si_;}

protected:

	SolverInterface * si_; /**< solver interface */

	double bestprimobj_; /**< best primal objective value */
	double bestdualobj_; /**< best dual objective value */

	int nworkers_;          /**< number of workers */
	int worker_;            /**< worker ID in communication */
	int nsubprobs_;         /**< number of subproblems for the current worker */
	int maxnumsubprobs_;    /**< maximum number of subproblems among workers */
	int * subindex_;        /**< array of subproblem indices */
	double * subprimobj_;   /**< subproblem primal objective values */
	double * subdualobj_;   /**< subproblem dual objective values */
	double ** subsolution_; /**< subproblem solution */
};

#endif /* SRC_SOLVER_DUALDECOMP_DDMASTER_H_ */
