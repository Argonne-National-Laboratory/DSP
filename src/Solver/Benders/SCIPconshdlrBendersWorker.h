/*
 * SCIPconshdlrBendersWorker.h
 *
 *  Created on: Feb 18, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVERINTERFACE_SCIPCONSHDLRBENDERSWORKER_H_
#define SRC_SOLVERINTERFACE_SCIPCONSHDLRBENDERSWORKER_H_

#include "Utility/DspMpi.h"
#include "Solver/Benders/SCIPconshdlrBenders.h"
#include "Solver/Benders/SCIPconshdlrBaseBendersWorker.h"

/** A class for implementing parallel Benders constraint handler */
class SCIPconshdlrBendersWorker : public SCIPconshdlrBenders, public SCIPconshdlrBaseBendersWorker
{
public:

	/** default constructor */
	SCIPconshdlrBendersWorker(SCIP *scip, int sepapriority, MPI_Comm comm)
		: SCIPconshdlrBenders(scip, "Benders", sepapriority),
		  SCIPconshdlrBaseBendersWorker(comm) {}

	/** default destructor */
	virtual ~SCIPconshdlrBendersWorker() {}

	/** set model pointer */
	virtual void setDecModel(DecModel *model)
	{
		model_ = model;
		initialize(model_->getNumSubproblems());
	}

protected:

	/** generate Benders cuts */
	virtual void generateCuts(
		int size,  /**< [in] size of x */
		double *x, /**< [in] master solution */
		OsiCuts *cuts /**< [out] cuts generated */)
	{
		generateCutsBase(model_->getNumSubproblems(), nvars_, naux_, x, probability_, cuts);
	}
};

#endif /* SRC_SOLVERINTERFACE_SCIPCONSHDLRBENDERSWORKER_H_ */
