/*
 * SolverInterfaceClp.h
 *
 *  Created on: Dec 9, 2014
 *      Author: kibaekkim
 */

#ifndef SOLVERINTERFACECLP_H_
#define SOLVERINTERFACECLP_H_

#include "SolverInterface/SolverInterfaceOsi.h"

class SolverInterfaceClp: public SolverInterfaceOsi
{
public:

	/** default constructor */
	SolverInterfaceClp(DspParams * par);

	/** copy constructor */
	SolverInterfaceClp(SolverInterfaceClp * si);

	/** copy constructor */
	SolverInterfaceClp(DspParams * par, OsiSolverInterface * si);

	/** default destructor */
	virtual ~SolverInterfaceClp() {}

protected:

	/** initialize solver interface */
	virtual STO_RTN_CODE initialize();
};

#endif /* SOLVERINTERFACECLP_H_ */
