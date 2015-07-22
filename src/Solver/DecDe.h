/*
 * DecDe.h
 *
 *  Created on: Sep 24, 2014
 *      Author: kibaekkim, ctjandra
 */

#ifndef DECDE_H_
#define DECDE_H_

#include <typeinfo>

/** DSP */
#include "Solver/DecSolver.h"
#include "Solver/SolverInterface.h"

class DecDe : public DecSolver {
public:

	/** default constructor */
	DecDe();

	/** default destructor */
	virtual ~DecDe();

	/** solve */
	virtual STO_RTN_CODE solve();

private:

	SolverInterface * si_; /**< my solver interface */
};

#endif /* DECDE_H_ */
