/*
 * BdMaster.h
 *
 *  Created on: Feb 18, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_BENDERS_BDMASTER_H_
#define SRC_SOLVER_BENDERS_BDMASTER_H_

/** DSP */
#include "Solver/DecSolver.h"
#include "Solver/Benders/SCIPconshdlrBenders.h"
#include "SolverInterface/SolverInterface.h"

class BdMaster: public DecSolver {

public:

	/** constructor */
	BdMaster(DspParams * par, DecModel * model, DspMessage * message);

	/** destructor */
	virtual ~BdMaster();

	/** initialize */
	virtual DSP_RTN_CODE init();

	/** solve */
	virtual DSP_RTN_CODE solve();

public:

	/** get solver interface pointer */
	virtual SolverInterface * getSiPtr() {return si_;}

public:

	/** set objective bounds */
	virtual DSP_RTN_CODE setObjectiveBounds(double upper, double lower);

	/** set constraint handler for Benders cut generation */
	virtual DSP_RTN_CODE setConshdlr(SCIPconshdlrBenders * conshdlr);

	/** set auxiliary column data */
	virtual DSP_RTN_CODE setAuxVarData(int size, double * obj, double * clbd, double * cubd);

	/** set initial solutions */
	virtual DSP_RTN_CODE setSolutions(Solutions initsols);

	/** set branching priorities */
	virtual DSP_RTN_CODE setBranchingPriority(
			int   size,      /**< size of array */
			int * priorities /**< branch priority */);

protected:

	/** create problem */
	virtual DSP_RTN_CODE createProblem();

protected:

	SolverInterface * si_; /**< my solver interface */

	int      naux_;     /**< number of auxiliary variables */
	double * obj_aux_;  /**< auxiliary variable objectives */
	double * clbd_aux_; /**< auxiliary variable lower bounds */
	double * cubd_aux_; /**< auxiliary variable upper bounds */
};

#endif /* SRC_SOLVER_BENDERS_BDMASTER_H_ */
