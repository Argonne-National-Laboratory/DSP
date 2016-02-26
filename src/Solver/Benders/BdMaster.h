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

	typedef vector<CoinPackedVector*> Solutions;

public:

	/** constructor */
	BdMaster(DspParams * par, DecModel * model, StoMessage * message);

	/** destructor */
	virtual ~BdMaster();

	/** initialize */
	virtual STO_RTN_CODE init();

	/** solve */
	virtual STO_RTN_CODE solve();

public:

	/** get solver interface pointer */
	virtual SolverInterface * getSiPtr() {return si_;}

public:

	/** set dual objective value */
	virtual STO_RTN_CODE setDualObjective(double dualobj);

	/** set constraint handler for Benders cut generation */
	virtual STO_RTN_CODE setConshdlr(SCIPconshdlrBenders * conshdlr);

	/** set auxiliary column data */
	virtual STO_RTN_CODE setAuxVarData(int size, double * obj, double * clbd, double * cubd);

	/** set initial solutions */
	virtual STO_RTN_CODE setSolutions(Solutions initsols);

	/** set branching priorities */
	virtual STO_RTN_CODE setBranchingPriority(
			int   size,      /**< size of array */
			int * priorities /**< branch priority */);

protected:

	/** create problem */
	virtual STO_RTN_CODE createProblem();

protected:

	SolverInterface * si_; /**< my solver interface */

	int      naux_;     /**< number of auxiliary variables */
	double * obj_aux_;  /**< auxiliary variable objectives */
	double * clbd_aux_; /**< auxiliary variable lower bounds */
	double * cubd_aux_; /**< auxiliary variable upper bounds */
};

#endif /* SRC_SOLVER_BENDERS_BDMASTER_H_ */
