/*
 * DwSolverSerial.h
 *
 *  Created on: Dec 5, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DANTZIGWOLFE_DWSOLVERSERIAL_H_
#define SRC_SOLVER_DANTZIGWOLFE_DWSOLVERSERIAL_H_

#include "Solver/DecSolver.h"
#include "TreeSearch/DspNodeSolution.h"
#include "Solver/DantzigWolfe/DwMaster.h"
#include "Solver/DantzigWolfe/DwModel.h"

class DwSolverSerial: public DecSolver {
public:

    /** default constructor */
	DwSolverSerial(
			DecModel *   model,  /**< model pointer */
			DspParams *  par,    /**< parameters */
			DspMessage * message /**< message pointer */);

	virtual DwSolverSerial* clone() const {return NULL;}

	virtual ~DwSolverSerial();

	/** initialize */
	virtual DSP_RTN_CODE init();

	/** solve */
    virtual DSP_RTN_CODE run() {return solve();}
	virtual DSP_RTN_CODE solve();

	/** finalize */
	virtual DSP_RTN_CODE finalize();

protected:

	DwMaster* master_;
	DwWorker* worker_;
	DwModel * alps_; /**< Alps model pointer */

};

#endif /* SRC_SOLVER_DANTZIGWOLFE_DWSOLVERSERIAL_H_ */
