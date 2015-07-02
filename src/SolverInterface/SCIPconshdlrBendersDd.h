/*
 * SCIPconshdlrBendersDd.h
 *
 *  Created on: Dec 14, 2014
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_SCIPCONSHDLRBENDERSDD_H_
#define SRC_SOLVER_SCIPCONSHDLRBENDERSDD_H_

/** SCIP */
#include "objscip/objconshdlr.h"

/** COIN */
#include "OsiCuts.hpp"

/** DSP */
#include "Solver/TssBdSub.h"

class SCIPconshdlrBendersDd : public scip::ObjConshdlr
{
public:

	enum whereFrom
	{
		from_scip_sepalp = 0,
		from_scip_sepasol,
		from_scip_enfolp,
		from_scip_enfops,
		from_scip_check
	};

	/** default constructor */
	SCIPconshdlrBendersDd(SCIP * scip) :
		ObjConshdlr(scip, "BendersDd", "Benders cuts",
				-2000000, /**< priority of the constraint handler for separation */
				-2000000, /**< priority of the constraint handler for constraint enforcing */
				-2000000, /**< priority of the constraint handler for checking infeasibility (and propagation) */
				1,        /**< always call separator */
				-1,       /**< disable constraint propagation */
				1,        /**< always use all constraints (no aging) */
				0,        /**< disable the preprocessing callback of the constraint handler */
				TRUE,     /**< delay separation method */
				FALSE,    /**< do not delay propatation method */
				FALSE,    /**< do not delay presolving method */
				TRUE,     /**< do not skip constraint handler even if no constraints are available. */
				SCIP_PROPTIMING_BEFORELP /**< propagation method is called before solving LP */),
		tss_(NULL), nvars_(0), vars_(NULL), nCutsReused_(0), ncols_(0), ncols_first_(0), obj_first_(NULL)
	{
		cutsToAdd_ = new OsiCuts;
		cutsAdded_ = new OsiCuts;
	}

	/** default constructor */
	virtual ~SCIPconshdlrBendersDd()
	{
		/** nothing to do */
	}

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

public:

	/** assign necessary data */
	void assignData(TssBdSub * tss, int ncols, int ncols_first, double * obj_first)
	{
		tss_ = tss;
		ncols_ = ncols;
		ncols_first_ = ncols_first;
		obj_first_ = obj_first;
	}

	/** set original variable pointers */
	SCIP_RETCODE setOriginalVariables(int nvars, SCIP_Var ** vars);

	/** get cutsAdded_ */
	const OsiCuts * getCutsAdded() {return cutsAdded_;}

	/** get cutsToAdd_ */
	const OsiCuts * getCutsToAdd() {return cutsToAdd_;}

	/** get number of cuts reused */
	int getNumCutsReused() {return nCutsReused_;}

	/** set cutsToAdd_ */
	void setCutsToAdd(OsiCuts * cuts);

	/** clean cutsAdded */
	void clearCutsAdded();

private:

	SCIP_RETCODE sepaBenders(
			SCIP * scip,
			SCIP_CONSHDLR * conshdlr,
			SCIP_SOL * sol,
			whereFrom where,
			SCIP_RESULT * result);

	/** construct upper bounding cuts */
	void constructCuts(OsiCuts & cs);

private:

	/** TODO: These could go in SCIP_ConsData. Does it matter? */
	TssBdSub * tss_;   /**< pointer to cut generator */
	int nvars_;        /**< number of original variables */
	SCIP_Var ** vars_; /**< pointer array to original variables */
	OsiCuts * cutsToAdd_; /**< cut pool */
	OsiCuts * cutsAdded_; /**< cut pool */
	int nCutsReused_; /**< number of cuts reused */

	int ncols_;          /**< number of Dd subproblem columns */
	int ncols_first_;    /**< number of first-stage columns */
	double * obj_first_; /**< first-stage objective coefficients */

};

/** creates and captures a Benders constraint */
SCIP_RETCODE SCIPcreateConsBendersDd(
		SCIP * scip,
		SCIP_CONS ** cons,
		const char * name);

#endif /* SRC_SOLVER_SCIPCONSHDLRBENDERSDD_H_ */
