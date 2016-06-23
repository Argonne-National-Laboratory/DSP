/*
 * SCIPconshdlrBendersMPI.h
 *
 *  Created on: Nov 20, 2015
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVERINTERFACE_SCIPCONSHDLRBENDERSMPI_H_
#define SRC_SOLVERINTERFACE_SCIPCONSHDLRBENDERSMPI_H_


/** MPI */
#include "mpi.h"

/** DSP */
#include "SCIPconshdlrBenders.h"

class SCIPconshdlrBendersMPI: public SCIPconshdlrBenders {
public:

	/** default constructor */
	SCIPconshdlrBendersMPI(SCIP * scip, int sepapriority, MPI_Comm comm):
		SCIPconshdlrBenders(scip,sepapriority),
		comm_(comm), comm_rank_(0), comm_size_(0),
		procIdxSize_(0), nSubs_(0), probability_(NULL),
		recvcounts_(NULL), displs_(NULL),
		cut_indices_(NULL), cut_status_(NULL), my_status_(NULL),
		cutval_(NULL), cutrhs_(NULL) {}

	/** default constructor */
	virtual ~SCIPconshdlrBendersMPI()
	{
		/** nothing to do */
	}

	/** initialize MPI related things */
	void initializeMpi(
			int n, /**< number of scenarios */
			int procIdxSize,
			const double * probability);

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

protected:

	SCIP_RETCODE sepaBenders(
			SCIP * scip,
			SCIP_CONSHDLR * conshdlr,
			SCIP_SOL * sol,
			whereFrom where,
			SCIP_RESULT * result);

private:

	/** aggregate cuts in cuts_in to cuts_out */
	void aggregateCuts(OsiCuts cuts_in, OsiCuts &cuts_out);

private:

	MPI_Comm comm_;
	int comm_rank_;
	int comm_size_;
	int procIdxSize_;
	int nSubs_;
	const double * probability_;

	/** cut communication */
	int * recvcounts_;
	int * displs_;

	int * cut_indices_; /**< scenario indices of cuts from MPI_Gather */
	int * cut_status_;  /**< cut generation status from MPI_Gatherv */
	int * my_status_;   /**< cut generation status of the current process */

	double ** cutval_;
	double * cutrhs_;
};

#endif /* SRC_SOLVERINTERFACE_SCIPCONSHDLRBENDERSMPI_H_ */
