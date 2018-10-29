/*
 * SolverInterfaceCpx.h
 *
 *  Created on: Mar 14, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVERINTERFACE_SOLVERINTERFACECPX_H_
#define SRC_SOLVERINTERFACE_SOLVERINTERFACECPX_H_

#include "cplex.h"
#include "SolverInterface/SolverInterfaceOsi.h"

/** Coin-OR */
#include "OsiCpxSolverInterface.hpp"

class SolverInterfaceCpx: public SolverInterfaceOsi {
public:

	/** constructor */
	SolverInterfaceCpx(DspParams * par);

	/** copy constructor */
	SolverInterfaceCpx(SolverInterfaceCpx * si);

	/** copy constructor */
	SolverInterfaceCpx(DspParams * par, OsiSolverInterface * si);

	/** default destructor */
	virtual ~SolverInterfaceCpx();

	/** solve */
	virtual void solve();

public:

	/** solution status */
	virtual DSP_RTN_CODE getStatus();

	/** get dual bound (lower bound in minimization) */
	virtual double getDualBound();

	/** set node limit */
	virtual void setNodeLimit(int limit);

	/** set iteration limit */
	virtual void setIterLimit(int limit);

	/** set wall time limit */
	virtual void setTimeLimit(double sec);

	/** set gap tolerance */
	virtual void setGapTol(double tol);

protected:

	/** initialize solver interface */
	virtual DSP_RTN_CODE initialize();

private:

	OsiCpxSolverInterface * cpx_;

public:

	bool useBarrier_; /**< use barrier solver */
};

#endif /* SRC_SOLVERINTERFACE_SOLVERINTERFACECPX_H_ */
