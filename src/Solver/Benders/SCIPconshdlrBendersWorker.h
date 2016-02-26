/*
 * SCIPconshdlrBendersWorker.h
 *
 *  Created on: Feb 18, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVERINTERFACE_SCIPCONSHDLRBENDERSWORKER_H_
#define SRC_SOLVERINTERFACE_SCIPCONSHDLRBENDERSWORKER_H_

/** MPI */
#include "mpi.h"

/** DSP */
#include "Solver/Benders/SCIPconshdlrBenders.h"

class SCIPconshdlrBendersWorker: public SCIPconshdlrBenders {
public:

	/** default constructor */
	SCIPconshdlrBendersWorker(SCIP * scip, int sepapriority, MPI_Comm comm);

	/** default destructor */
	virtual ~SCIPconshdlrBendersWorker();

	/** destructor of constraint handler to free user data (called when SCIP is exiting) */
	virtual SCIP_DECL_CONSFREE(scip_free);

	/** clone method which will be used to copy constraint handler and variable pricer objects */
	virtual SCIP_DECL_CONSHDLRCLONE(ObjProbCloneable* clone);

	/** set number of subproblems */
	virtual void setNumSubprobs(int nsubprobs);

protected:

	/** generate Benders cuts */
	virtual void generateCuts(
			int size,      /**< [in] size of x */
			double * x,    /**< [in] master solution */
			int where,     /**< [in] where to be called */
			OsiCuts * cuts /**< [out] cuts generated */);

private:

	MPI_Comm comm_;
	int comm_rank_;
	int comm_size_;

	/** cut communication */
	int * recvcounts_;
	int * displs_;

	int * cut_indices_;          /**< scenario indices of cuts from MPI_Gather */
	int * cut_status_;           /**< cut generation status from MPI_Gatherv */
};

#endif /* SRC_SOLVERINTERFACE_SCIPCONSHDLRBENDERSWORKER_H_ */
