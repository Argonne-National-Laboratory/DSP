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
	SCIPconshdlrIntBendersWorker(SCIP *scip, int sepapriority, int param_seps_solver, MPI_Comm comm)
		: SCIPconshdlrIntBenders(scip, "Benders", sepapriority, param_seps_solver),
		  SCIPconshdlrBaseBendersWorker(comm) {}

	/** default destructor */
	virtual ~SCIPconshdlrIntBendersWorker() {}

	/** set model pointer */
	virtual void setDecModel(DecModel *model);

protected:
	/** generate Benders cuts */
	virtual void generateCuts(
		int size,  /**< [in] size of x */
		double *x, /**< [in] master solution */
		OsiCuts *cuts /**< [out] cuts generated */)
	{
		generateCutsBase(model_->getNumSubproblems(), nvars_, naux_, x, probability_, cuts);
	}

	/** evaluate recourse */
	virtual SCIP_RETCODE evaluateRecourse(
		SCIP *scip,	   /**< [in] scip pointer */
		SCIP_SOL *sol, /**< [in] solution to evaluate */
		double *values /**< [out] evaluated recourse values */);
};

#endif /* SRC_SOLVERINTERFACE_SCIPCONSHDLRINTBENDERSWORKER_H_ */
