/*
 * DdMaster.h
 *
 *  Created on: Feb 9, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DUALDECOMP_DDMASTER_H_
#define SRC_SOLVER_DUALDECOMP_DDMASTER_H_

#include "Solver/DecSolver.h"

/** Base class for the decomposition master solver */
class DdMaster: public DecSolver {

	friend class DdMW;
	friend class DdMWSerial;
	friend class DdMWAsync;
	friend class DdMWSync;

public:

	/** A default constructor. */
	DdMaster(
			DecModel *   model,  /**< model pointer */
			DspParams *  par,    /**< parameter pointer */
			DspMessage * message /**< message pointer */);

	/** A copy constructor */
	DdMaster(const DdMaster& rhs);

	/** A default destructor. */
	virtual ~DdMaster();

	/** A clone function */
	virtual DdMaster* clone() const {
		return new DdMaster(*this);
	}

	/** A virtual member for initializing solver. */
	virtual DSP_RTN_CODE init();

	/** A virtual member for solving problem. */
	virtual DSP_RTN_CODE solve() {return DSP_RTN_OK;}

	/** A virtual memeber for finalizing solver. */
	virtual DSP_RTN_CODE finalize() {return DSP_RTN_OK;}

	/** A virtual member for updating problem */
	virtual DSP_RTN_CODE updateProblem() {return DSP_RTN_OK;}

	/** A virtual member for setting initial solution */
	virtual DSP_RTN_CODE setInitSolution(const double * sol);

	/** A virtual member for termination test */
	virtual DSP_RTN_CODE terminationTest() {return status_;}

	/** A const member to return lambda */
	const double * getLambda() {return lambda_;}

protected:

	/** create DspOsi for master */
	virtual DspOsi * createDspOsi();

	const double* lambda_; /**< pointer to the lambda part (Lagrangian multiplier with respect to the relaxed constraints) of the solution */
	std::vector<double> subprimobj_;   /**< subproblem primal objective values */
	std::vector<double> subdualobj_;   /**< subproblem dual objective values */
	double ** subsolution_; /**< subproblem solution */
};

#endif /* SRC_SOLVER_DUALDECOMP_DDMASTER_H_ */
