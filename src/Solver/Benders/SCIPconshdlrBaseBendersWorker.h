/*
 * SCIPconshdlrBaseBendersWorker.h
 */

#ifndef SRC_SOLVERINTERFACE_SCIPCONSHDLRBASEBENDERSWORKER_H_
#define SRC_SOLVERINTERFACE_SCIPCONSHDLRBASEBENDERSWORKER_H_

#include "Utility/DspMpi.h"

/** A base class for implementing parallel Benders constraint handler */
class SCIPconshdlrBaseBendersWorker
{
public:
	/** default constructor */
	SCIPconshdlrBaseBendersWorker(MPI_Comm comm);

	/** default destructor */
	virtual ~SCIPconshdlrBaseBendersWorker();

	/** initialize */
	virtual void initialize(int nsubprobs);

protected:
	/** generate Benders cuts */
	virtual void generateCutsBase(
		int nsubprobs, /**< [in] number of subproblems */
		int nvars,	   /**< [in] number of variables */
		int naux,	   /**< [in] number of auxiliary variables */
		double *x,	   /**< [in] master solution */
		OsiCuts *cuts /**< [out] cuts generated */);

	/** generate Benders cuts */
	virtual void aggregateCutsBase(
		int nsubprobs,	 /**< [in] number of subproblems */
		int nvars,		 /**< [in] number of variables */
		int naux,		 /**< [in] number of auxiliary variables */
		double **cutvec, /**< [in] cut vector */
		double *cutrhs,	 /**< [in] cut right-hand side */
		OsiCuts *cuts /**< [out] cuts generated */);

	MPI_Comm comm_;
	int comm_rank_;
	int comm_size_;

	/** cut communication */
	int *recvcounts_;
	int *displs_;
	int *cut_index_;  /**< subproblem index for which cut is generated */
	int *cut_status_; /**< cut generation status from MPI_Gatherv */
};

#endif /* SRC_SOLVERINTERFACE_SCIPCONSHDLRBASEBENDERSWORKER_H_ */
