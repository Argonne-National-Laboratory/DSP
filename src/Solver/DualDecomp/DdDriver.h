/*
 * DdDriver.h
 *
 *  Created on: Feb 16, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DUALDECOMP_DDDRIVER_H_
#define SRC_SOLVER_DUALDECOMP_DDDRIVER_H_

#include "Solver/DecSolver.h"
#include "Solver/DualDecomp/DdMW.h"

/** A base driver class for dual decomposition */
class DdDriver: public DecSolver {
public:

	/** A default constructor. */
	DdDriver(
			DecModel *   model,  /**< model pointer */
			DspParams *  par,    /**< parameters */
			DspMessage * message /**< message pointer */);

	/** A copy constructor. */
	DdDriver(const DdDriver& rhs);

	/** A default destructor. */
	virtual ~DdDriver();

	/** A clone function */
	virtual DdDriver* clone() const {
		return new DdDriver(*this);
	}

	/** A virtual member for initializing the driver. */
	virtual DSP_RTN_CODE init() {return DSP_RTN_OK;}

	/** A virtual member for running the driver. */
	virtual DSP_RTN_CODE run() {return DSP_RTN_OK;}
	virtual DSP_RTN_CODE solve() {return DSP_RTN_ERR;}

	/** A virtual memeber for finalizing solver. */
	virtual DSP_RTN_CODE finalize() {return DSP_RTN_OK;}

protected:

	DdMW * mw_;
};

#endif /* SRC_SOLVER_DUALDECOMP_DDDRIVER_H_ */
