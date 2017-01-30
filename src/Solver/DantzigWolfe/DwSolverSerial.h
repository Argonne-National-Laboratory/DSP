/*
 * DwSolverSerial.h
 *
 *  Created on: Dec 5, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DANTZIGWOLFE_DWSOLVERSERIAL_H_
#define SRC_SOLVER_DANTZIGWOLFE_DWSOLVERSERIAL_H_

#include <DecSolver.h>
#include <DantzigWolfe/DwMaster.h>
#include "TreeSearch/DspModel.h"

class DwSolverSerial: public DecSolver {
public:

    /** default constructor */
	DwSolverSerial(
			DecModel *   model,  /**< model pointer */
			DspParams *  par,    /**< parameters */
			DspMessage * message /**< message pointer */);

	virtual ~DwSolverSerial();

	/** initialize */
	virtual DSP_RTN_CODE init();

	/** solve */
	virtual DSP_RTN_CODE solve();

	/** finalize */
	virtual DSP_RTN_CODE finalize();

protected:

	DwMaster* master_;
	DwWorker* worker_;
	DspModel * alps_; /**< Alps model pointer */

};

#endif /* SRC_SOLVER_DANTZIGWOLFE_DWSOLVERSERIAL_H_ */
