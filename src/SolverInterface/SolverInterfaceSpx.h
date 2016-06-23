/*
 * SolverInterfaceSpx.h
 *
 *  Created on: Feb 9, 2015
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_SOLVERINTERFACESPX_H_
#define SRC_SOLVER_SOLVERINTERFACESPX_H_

#include "SolverInterface/SolverInterfaceOsi.h"

class SolverInterfaceSpx: public SolverInterfaceOsi
{
public:

	/** default constructor */
	SolverInterfaceSpx(DspParams * par);

	/** copy constructor */
	SolverInterfaceSpx(SolverInterfaceSpx * si);

	/** copy constructor */
	SolverInterfaceSpx(DspParams * par, OsiSolverInterface * si);

	/** default destructor */
	virtual ~SolverInterfaceSpx() {}

protected:

	/** initialize solver interface */
	virtual DSP_RTN_CODE initialize();
};

#endif /* SRC_SOLVER_SOLVERINTERFACESPX_H_ */
