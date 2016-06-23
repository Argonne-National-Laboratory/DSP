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
#include "SolverInterface/SCIPconshdlrBenders.h"

class SCIPconshdlrBendersWorker: public SCIPconshdlrBenders {
public:

	/** default constructor */
	SCIPconshdlrBendersWorker(SCIP * scip, int sepapriority, MPI_Comm comm):
		comm_(comm), recvcounts_(NULL), displs_(NULL), probability_(NULL),
		cut_indices_(NULL), cut_status_(NULL)
	{
		MPI_Comm_size(comm_, &comm_size_);
		MPI_Comm_rank(comm_, &comm_rank_);
	}

	/** default destructor */
	virtual ~SCIPconshdlrBendersWorker();

	/** destructor of constraint handler to free user data (called when SCIP is exiting) */
	virtual SCIP_DECL_CONSFREE(scip_free);

	/** transforms constraint data into data belonging to the transformed problem */
	virtual SCIP_DECL_CONSTRANS(scip_trans);

	/** separation method of constraint handler for LP solution */
	virtual SCIP_DECL_CONSSEPALP(scip_sepalp);

	/** separation method of constraint handler for arbitrary primal solution */
	virtual SCIP_DECL_CONSSEPASOL(scip_sepasol);

	/** constraint enforcing method of constraint handler for LP solutions */
	virtual SCIP_DECL_CONSENFOLP(scip_enfolp);

	/** constraint enforcing method of constraint handler for pseudo solutions */
	virtual SCIP_DECL_CONSENFOPS(scip_enfops);

	/** feasibility check method of constraint handler for primal solutions */
	virtual SCIP_DECL_CONSCHECK(scip_check);

	/** variable rounding lock method of constraint handler */
	virtual SCIP_DECL_CONSLOCK(scip_lock);

	/** returns whether the objective plugin is copyable */
	virtual SCIP_DECL_CONSHDLRISCLONEABLE(iscloneable) {return true;}

	/** clone method which will be used to copy constraint handler and variable pricer objects */
	virtual SCIP_DECL_CONSHDLRCLONE(ObjProbCloneable* clone);

	/** initialize */
	void init(
			int nsubprobs,                    /**< number of subproblems */
			const double * probability = NULL /**< probability of each subproblem */);

protected:

	SCIP_RETCODE sepaBenders(
			SCIP * scip,
			SCIP_CONSHDLR * conshdlr,
			SCIP_SOL * sol,
			whereFrom where,
			SCIP_RESULT * result);

	/** aggregate cuts in cuts_in to cuts_out */
	void aggregateCuts(OsiCuts cuts_in, OsiCuts &cuts_out);

private:

	MPI_Comm comm_;
	int comm_rank_;
	int comm_size_;

	/** cut communication */
	int * recvcounts_;
	int * displs_;

	const double * probability_; /**< probability of each subproblem */
	int * cut_indices_;          /**< scenario indices of cuts from MPI_Gather */
	int * cut_status_;           /**< cut generation status from MPI_Gatherv */
};

#endif /* SRC_SOLVERINTERFACE_SCIPCONSHDLRBENDERSWORKER_H_ */
