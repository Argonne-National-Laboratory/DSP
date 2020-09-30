/*
 * SCIPconshdlrIntBendersWorker.h
 */

#ifndef SRC_SOLVERINTERFACE_SCIPCONSHDLRINTBENDERSWORKER_H_
#define SRC_SOLVERINTERFACE_SCIPCONSHDLRINTBENDERSWORKER_H_

#include "Solver/Benders/SCIPconshdlrIntBenders.h"
#include "Solver/Benders/SCIPconshdlrBaseBendersWorker.h"

/** A class for implementing parallel integer Benders constraint handler */
class SCIPconshdlrIntBendersWorker : public SCIPconshdlrIntBenders, public SCIPconshdlrBaseBendersWorker
{
public:
	/** default constructor */
	SCIPconshdlrIntBendersWorker(SCIP *scip, int sepapriority, MPI_Comm comm)
		: SCIPconshdlrIntBenders(scip, "Benders", sepapriority),
		  SCIPconshdlrBaseBendersWorker(comm) {}

	/** default destructor */
	virtual ~SCIPconshdlrIntBendersWorker() {}

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
		generateCutsBase(model_->getNumSubproblems(), nvars_, naux_, x, cuts);
	}

	/** evaluate recourse: 
	 * This function is called only when the recourse has integer variables.
	*/
	virtual SCIP_RETCODE evaluateRecourse(
		SCIP *scip,	   /**< [in] scip pointer */
		SCIP_SOL *sol, /**< [in] solution to evaluate */
		double *values /**< [out] evaluated recourse values */);
};

#endif /* SRC_SOLVERINTERFACE_SCIPCONSHDLRINTBENDERSWORKER_H_ */
