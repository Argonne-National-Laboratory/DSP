/*
 * DdSolver.h
 *
 *  Created on: Feb 9, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DUALDECOMP_DDSOLVER_H_
#define SRC_SOLVER_DUALDECOMP_DDSOLVER_H_

#include "Solver/DecSolver.h"

class DdSolver: public DecSolver {
public:

	/** constructor */
	DdSolver(DspParams * par, DecModel * model, DspMessage * message):
		DecSolver(par, model, message) {}

	/** destructor */
	virtual ~DdSolver() {}

	/** solve */
	virtual DSP_RTN_CODE solve() {return DSP_RTN_OK;}

protected:

	/** create problem */
	virtual DSP_RTN_CODE createProblem() {return DSP_RTN_OK;}
};

#endif /* SRC_SOLVER_DUALDECOMP_DDSOLVER_H_ */
