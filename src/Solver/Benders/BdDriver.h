/*
 * BdDriver.h
 *
 *  Created on: Feb 17, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_BENDERS_BDDRIVER_H_
#define SRC_SOLVER_BENDERS_BDDRIVER_H_

#include "Solver/DspDriver.h"
#include "Solver/Benders/BdMW.h"
#include "Solver/Benders/BdMaster.h"
#include "Solver/Benders/BdWorker.h"
#include "Solver/Benders/SCIPconshdlrBenders.h"

class BdDriver: public DspDriver {

	typedef vector<CoinPackedVector*> Solutions;

public:

	/** constructor */
	BdDriver(DspParams * par, DecModel * model);

	/** constructor */
	BdDriver(DspParams * par, DecModel * model, MPI_Comm comm);

	/** destructor */
	virtual ~BdDriver();

	/** initialize */
	virtual STO_RTN_CODE init();

	/** run */
	virtual STO_RTN_CODE run();

public:

	/** set dual objective */
	virtual void setDualObjective(double dualobj) {dualobj_ = dualobj;}

	/** set auxiliary variable data */
	virtual STO_RTN_CODE setAuxVarData(int size, double* obj, double* clbd, double* cubd);

	/** set initial solution */
	virtual STO_RTN_CODE setSolution(
			int      size,    /**< size of array */
			double * solution /**< solution */);

	/** set branch priorities */
	virtual STO_RTN_CODE setPriorities(
			int   size,      /**< size of array */
			int * priorities /**< branch priority */);

private:

	/** find lower bound */
	STO_RTN_CODE findLowerBound();

	/** constraint handler */
	SCIPconshdlrBenders * constraintHandler();

private:

	bool useMPI_; /**< indicating whether to use MPI library or not */

	MPI_Comm comm_;
	int comm_rank_;
	int comm_size_;

	BdMW * mw_;
	BdMaster * master_;
	BdWorker * worker_;

private:

	Solutions initsols_; /**< set of initial solutions */
	int numPriorities_;  /**< length of branch priorities */
	int * priorities_;   /**< branch priority */
};

#endif /* SRC_SOLVER_BENDERS_BDDRIVER_H_ */
