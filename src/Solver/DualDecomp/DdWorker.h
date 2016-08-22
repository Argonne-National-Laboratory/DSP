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

	enum {
		Base = 0,
		LB,  /**< lower bounder */
		UB,  /**< upper bounder */
		CGBd /**< Benders cut generator */
	};

	/** constructor */
	DdWorker(DspParams * par, DecModel * model, DspMessage * message);

	/** destructor */
	virtual ~DdWorker();

	virtual int getType() {return Base;}

};

#endif /* SRC_SOLVER_DUALDECOMP_DDWORKER_H_ */
