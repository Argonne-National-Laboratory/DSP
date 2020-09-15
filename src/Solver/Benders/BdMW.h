/*
 * BdMW.h
 *
 *  Created on: Feb 17, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_BENDERS_BDMW_H_
#define SRC_SOLVER_BENDERS_BDMW_H_

#include "Solver/BaseMasterWorker.h"
#include "Solver/Benders/BdMaster.h"
#include "Solver/Benders/BdWorker.h"
#include "Solver/Benders/SCIPconshdlrBenders.h"
#include "Solver/Benders/SCIPheurIntBenders.h"

/** A base class for the Benders master-worker framework */
class BdMW: public BaseMasterWorker {

public:

	enum
	{
		MASTER_NEEDS_CUTS = 0,
		MASTER_STOPPED
	};

	/** A default constructor. */
	BdMW(
			DecModel *   model,  /**< model pointer */
			DspParams *  par,    /**< parameters */
			DspMessage * message /**< message pointer */);

	/** A copy constructor. */
	BdMW(const BdMW& rhs);

	/** A default destructor. */
	virtual ~BdMW();

	/** A clone function */
	virtual BdMW* clone() const {
		return new BdMW(*this);
	}

	/** initialize */
	virtual DSP_RTN_CODE init() {return DSP_RTN_OK;}

	/** run the framework */
	virtual DSP_RTN_CODE run() {return DSP_RTN_OK;}

	/** finalize */
	virtual DSP_RTN_CODE finalize() {return DSP_RTN_OK;}

	virtual void warning_relaxation() {
		message_->print(0, "************************************************************************\n");
		message_->print(0, "* WARNING: Integer variables are detected but relaxed in the recourse. *\n");
		message_->print(0, "************************************************************************\n");
	}

protected:

	/** run master process */
	virtual DSP_RTN_CODE runMaster() {return DSP_RTN_OK;}

	/** run worker processes */
	virtual DSP_RTN_CODE runWorker() {return DSP_RTN_OK;}

	/** constraint handler */
	virtual SCIPconshdlrBenders * constraintHandler() {return NULL;}

	/** heuristic handler */
	virtual SCIPheurIntBenders * heuristicHandler() {return NULL;}

public:

	/** master pointer */
	BdMaster * getMasterPtr() {return master_;}

	/** Worker pointer */
	BdWorker * getWorkerPtr() {return worker_;}

protected:

	DecModel * model_;     /**< DecModel object */
	DspParams * par_;      /**< parameters */
	DspMessage * message_; /**< message */

	BdMaster * master_; /**< master */
	BdWorker * worker_; /**< worker */
};

#endif /* SRC_SOLVER_BENDERS_BDMW_H_ */
