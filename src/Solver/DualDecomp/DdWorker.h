/*
 * DdWorker.h
 *
 *  Created on: Feb 9, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DUALDECOMP_DDWORKER_H_
#define SRC_SOLVER_DUALDECOMP_DDWORKER_H_

#include "Solver/DecSolver.h"

/** An abstract class for decomposition worker */
class DdWorker : public DecSolver {
public:

	/** Type of workers */
	enum {
		Base = 0,
		LB,  /**< lower bounder */
		UB,  /**< upper bounder */
		CG,  /**< abstract cut generator */
		CGBd /**< Benders cut generator */
	};

	/** A default constructor. */
	DdWorker(
			DecModel *   model,  /**< model pointer */
			DspParams *  par,    /**< parameter pointer */
			DspMessage * message /**< message pointer */) :
	DecSolver(model, par, message) {}

	/** A copy constructor. */
	DdWorker(const DdWorker& rhs) : DecSolver(rhs) {}

	/** A default destructor. */
	virtual ~DdWorker() {}

	/** A clone function */
	virtual DdWorker* clone() const {
		return new DdWorker(*this);
	}

	/** A virtual member for initializing solver. */
	virtual DSP_RTN_CODE init() {return DSP_RTN_OK;}

	/** A virtual member for solving problem. */
	virtual DSP_RTN_CODE solve() {return DSP_RTN_OK;}

	/** A virtual memeber for finalizing solver. */
	virtual DSP_RTN_CODE finalize() {return DSP_RTN_OK;}

	/** A virtual member to return the worker type. */
	virtual int getType() {return Base;}
};

#endif /* SRC_SOLVER_DUALDECOMP_DDWORKER_H_ */
