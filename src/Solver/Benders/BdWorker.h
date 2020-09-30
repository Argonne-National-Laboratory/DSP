/*
 * BdWorker.h
 *
 *  Created on: Feb 18, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_BENDERS_BDWORKER_H_
#define SRC_SOLVER_BENDERS_BDWORKER_H_

#include "OsiCuts.hpp"
#include "Solver/Benders/BdSub.h"
#include "Model/DecModel.h"
#include "Utility/DspMessage.h"

/** A base class for Benders worker. This solves subproblem. */
class BdWorker {
public:

	/** A default constructor. */
	BdWorker(
		DecModel *   model,  /**< model pointer */
		DspParams *  par,    /**< parameters */
		DspMessage * message /**< message pointer */);

	/** A copy constructor. */
	BdWorker(const BdWorker& rhs);

	/** A default destructor. */
	virtual ~BdWorker();

	/** A clone function. */
	virtual BdWorker* clone() const {
		return new BdWorker(*this);
	}

	/** get BdSub pointer */
	virtual BdSub * getBdSubPtr() {return bdsub_;}

protected:

	DecModel * model_;     /**< DecModel object */
	DspParams * par_;      /**< parameters */
	DspMessage * message_; /**< message */

	BdSub * bdsub_; /**< benders subproblem */

	int parProcIdxSize_; /**< number of subproblems for this worker */
	int * parProcIdx_;   /**< subproblem indices for this worker */
};

#endif /* SRC_SOLVER_BENDERS_BDWORKER_H_ */
