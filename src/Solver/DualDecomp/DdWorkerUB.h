/*
 * DdWorkerUB.h
 *
 *  Created on: Mar 28, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DUALDECOMP_DDWORKERUB_H_
#define SRC_SOLVER_DUALDECOMP_DDWORKERUB_H_

#include <Solver/DualDecomp/DdWorkerLB.h>

class DdWorkerUB: public DdWorkerLB {
public:

	/** constructor */
	DdWorkerUB(DspParams * par, DecModel * model, DspMessage * message);

	/** destructor */
	virtual ~DdWorkerUB();

	/** solve */
	virtual DSP_RTN_CODE solve();

	virtual int getType() {return UB;}

	/** fix coupling variable values */
	DSP_RTN_CODE fixCouplingVariableValues(CoinPackedVector * val);
};

#endif /* SRC_SOLVER_DUALDECOMP_DDWORKERUB_H_ */