/*
 * SCIPconshdlrBenders.h
 *
 *  Created on: Dec 8, 2014
 *      Author: kibaekkim
 */

#ifndef SCIPCONSHDLRBENDERS_H_
#define SCIPCONSHDLRBENDERS_H_

#include "scip/def.h"
#include "objscip/objconshdlr.h"
#include "Model/DecModel.h"
#include "Solver/Benders/BdSub.h"
#include "OsiCuts.hpp"

/** A base class for implementing Benders constraint handler */
class SCIPconshdlrBenders : public scip::ObjConshdlr
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
	SCIPconshdlrBenders(SCIP * scip, const char * name, int sepapriority) :
		ObjConshdlr(scip, name, "Benders cuts",
				sepapriority, /**< priority of the constraint handler for separation */
				sepapriority, /**< priority of the constraint handler for constraint enforcing */
				sepapriority, /**< priority of the constraint handler for checking infeasibility (and propagation) */
				-1,       /**< never call separator */
				-1,       /**< disable constraint propagation */
				1,        /**< always use all constraints (no aging) */
				0,        /**< disable the preprocessing callback of the constraint handler */
				TRUE,     /**< delay separation method */
				FALSE,    /**< do not delay propatation method */
#if SCIP_VERSION < 320
				FALSE,    /**< do not delay presolving method */
				TRUE,     /**< skip constraint handler even if no constraints are available. */
				SCIP_PROPTIMING_BEFORELP /**< propagation method is called before solving LP */),
#else
				TRUE,     /**< skip constraint handler even if no constraints are available. */
				SCIP_PROPTIMING_BEFORELP, /**< propagation method is called before solving LP */
				SCIP_PRESOLTIMING_FAST),
#endif
	model_(NULL),
	bdsub_(NULL),
	nvars_(0),
	vars_(NULL),
	naux_(0),
	probability_(NULL)
	{
		/** nothing to do */
	}

	/** default constructor */
	virtual ~SCIPconshdlrBenders()
	{
		/** nothing to do */
	}

	/** destructor of constraint handler to free user data (called when SCIP is exiting) */
	virtual SCIP_DECL_CONSFREE(scip_free);

	/** transforms constraint data into data belonging to the transformed problem */
	virtual SCIP_DECL_CONSTRANS(scip_trans);

	/** constraint enforcing method of constraint handler for LP solutions */
	virtual SCIP_DECL_CONSENFOLP(scip_enfolp);

	/** constraint enforcing method of constraint handler for pseudo solutions */
	virtual SCIP_DECL_CONSENFOPS(scip_enfops);

	/** feasibility check method of constraint handler for primal solutions */
	virtual SCIP_DECL_CONSCHECK(scip_check);

	/** variable rounding lock method of constraint handler */
	virtual SCIP_DECL_CONSLOCK(scip_lock);

	/** returns whether the objective plugin is copyable */
	virtual SCIP_DECL_CONSHDLRISCLONEABLE(iscloneable) {return false;}

public:

	/** get pointer to cut generator */
	virtual const BdSub * getBdSub() {return bdsub_;}

	/** get number of original variables, including auxiliary variables */
	virtual int getNumVars() {return nvars_;}

	/** get pointer to original variables */
//	virtual const SCIP_Var ** getVars() {return vars_;}

	/** get number of auxiliary variables */
	virtual int getNumAuxVars() {return naux_;}

	virtual bool isStochastic() {
		if (model_ && model_->isStochastic()) return true;
		else return false;
	}

	virtual bool isDro() {
		if (model_ && model_->isDro()) return true;
		else return false;
	}

public:

	/** set model pointer */
	virtual void setDecModel(DecModel * model) {model_ = model;}

	/** set pointer to cut generator */
	virtual void setBdSub(BdSub * bdsub);

	/** set original variable pointers */
	virtual SCIP_RETCODE setOriginalVariables(
			int nvars,        /**< number of original variables, including auxiliary variables */
			SCIP_Var ** vars, /**< original variables, including auxiliary variables */
			int         naux  /**< number of auxiliary variables */);

protected:

	virtual SCIP_RETCODE sepaBenders(
			SCIP * scip,
			SCIP_CONSHDLR * conshdlr,
			SCIP_SOL * sol,
			whereFrom where,
			SCIP_RESULT * result);

	/** generate Benders cuts */
	virtual void generateCuts(
			int size,      /**< [in] size of x */
			double * x,    /**< [in] master solution */
			int where,     /**< [in] where to be called */
			OsiCuts * cuts /**< [out] cuts generated */);

	/** generate Benders cuts */
	virtual void aggregateCuts(
			double ** cutvec, /**< [in] cut vector */
			double *  cutrhs, /**< [in] cut right-hand side */
			OsiCuts * cuts    /**< [out] cuts generated */);

	/** Has integer variables in the recouse? */
	virtual bool isIntegralRecourse() {
		return (bdsub_ != NULL && bdsub_->has_integer());
	}

	/** evaluate recourse: 
	 * This function is called only when the recourse has integer variables.
	*/
	virtual SCIP_RETCODE evaluateRecourse(
			SCIP * scip,    /**< [in] scip pointer */
			SCIP_SOL * sol, /**< [in] solution to evaluate */
			double * values /**< [out] evaluated recourse values */);

	/** compuate probability (used for DRO) */
	virtual void computeProbability(
			const double* recourse /**< [in] recourse values */);

	virtual SCIP_RETCODE addNoGoodCut(
			SCIP * scip,             /**< [in] scip pointer */
   			SCIP_CONSHDLR* conshdlr, /**< [in] constraint handler that creates the row */
			SCIP_RESULT * result     /**< [out] result */);

	virtual SCIP_RETCODE addIntOptimalityCut(
			SCIP * scip,             /**< [in] scip pointer */
   			SCIP_CONSHDLR* conshdlr, /**< [in] constraint handler that creates the row */
			double exact_recourse,   /**< [in] exact recourse value */
			double recourse_lb,      /**< [in] lower bound of recourse value */
			SCIP_RESULT * result     /**< [out] result */);

protected:

	DecModel *  model_;            /**< DecModel object */
	BdSub *     bdsub_;            /**< pointer to cut generator */
	int         nvars_;            /**< number of original variables */
	SCIP_Var ** vars_;             /**< pointer array to original variables */
	int         naux_;             /**< number of auxiliary variables */
	double*     probability_;      /**< array of probability */
};

/** creates and captures a Benders constraint */
SCIP_RETCODE SCIPcreateConsBenders(
		SCIP * scip,
		SCIP_CONS ** cons,
		const char * name);

#endif /* SCIPCONSHDLRBENDERS_H_ */
