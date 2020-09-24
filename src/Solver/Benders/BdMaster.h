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
#include "Solver/Benders/BdWorker.h"

/** A class for implementing the Benders master solver */
class BdMaster: public DecSolver {

public:

	/** A default constructor. */
	BdMaster(
			DecModel *   model,   /**< model pointer */
			DspParams *  par,     /**< parameter pointer */
			DspMessage * message /**< message pointer */);

	/** A copy constructor. */
	BdMaster(const BdMaster& rhs);

	/** A default destructor. */
	virtual ~BdMaster();

	/** A clone function */
	virtual BdMaster* clone() const {
		return new BdMaster(*this);
	}

	/** A virtual member for initializing solver. */
	virtual DSP_RTN_CODE init();

	/** A virtual member for solving problem. */
	virtual DSP_RTN_CODE solve();

	/** A virtual memeber for finalizing solver. */
	virtual DSP_RTN_CODE finalize() {return DSP_RTN_OK;}

	/**@name Set functions */
	//@{

	/** set worker pointer */
	void setWorkerPtr(BdWorker* worker) {worker_ = worker;}

	/** set objective bounds */
	virtual DSP_RTN_CODE setObjectiveBounds(double upper, double lower);

	/** set constraint handler for Benders cut generation */
	virtual DSP_RTN_CODE setConshdlr(SCIPconshdlrBenders * conshdlr);

	/** set auxiliary column data */
	virtual DSP_RTN_CODE setAuxVarData(int size, double * obj, double * clbd, double * cubd);

	/** set initial solutions */
	virtual DSP_RTN_CODE setSolutions(Solutions initsols);

	//@}

	virtual int getNumAuxVars() {return naux_;}

	virtual bool is_binary() {return is_binary_;}

protected:

	/** create problem */
	virtual DSP_RTN_CODE createProblem();

protected:

	BdWorker* worker_; /**< Benders worker */

	int      naux_;     /**< number of auxiliary variables */
	double * obj_aux_;  /**< auxiliary variable objectives */
	double * clbd_aux_; /**< auxiliary variable lower bounds */
	double * cubd_aux_; /**< auxiliary variable upper bounds */
	bool is_binary_;    /**< indicate whether the first stage is a pure binary or not */
};

#endif /* SRC_SOLVER_BENDERS_BDMASTER_H_ */
