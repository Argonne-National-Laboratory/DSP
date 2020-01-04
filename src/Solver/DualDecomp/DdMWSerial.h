/*
 * DdMWSerial.h
 *
 *  Created on: Apr 8, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DUALDECOMP_DDMWSERIAL_H_
#define SRC_SOLVER_DUALDECOMP_DDMWSERIAL_H_

#include "Solver/DualDecomp/DdMW.h"

/** A master-worker class for serial dual decomposition */
class DdMWSerial: public DdMW {
public:

	/** A default constructor. */
	DdMWSerial(
			DecModel *   model,  /**< model pointer */
			DspParams *  par,    /**< parameters */
			DspMessage * message /**< message pointer */);

	/** A copy constructor. */
	DdMWSerial(const DdMWSerial& rhs);

	/** A default destructor. */
	virtual ~DdMWSerial();

	/** A clone function */
	virtual DdMWSerial* clone() const {
		return new DdMWSerial(*this);
	}

	/** A virtual member for initializing the framework. */
	virtual DSP_RTN_CODE init();

	/** A virtual member for running the framework. */
	virtual DSP_RTN_CODE run();

	/** A virtual memeber for finalizing the framework. */
	virtual DSP_RTN_CODE finalize();

protected:

#ifdef DSP_HAS_SCIP
	/** generate Benders cuts */
	DSP_RTN_CODE generateBendersCuts(
			DdWorkerCGBd * workercg, /**< CG worker pointer */
			Solutions solutions, /**< solutions at which cuts are generated */
			OsiCuts & cuts       /**< cuts generated */);
#endif

};

#endif /* SRC_SOLVER_DUALDECOMP_DDMWSERIAL_H_ */
