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
	/** default constructor */
	SCIPconshdlrBenders(SCIP *scip, const char *name, int sepapriority);

	/** default constructor */
	virtual ~SCIPconshdlrBenders();

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

	/** separation method of constraint handler for LP solution */
	virtual SCIP_DECL_CONSSEPALP(scip_sepalp);

	/** separation method of constraint handler for arbitrary primal solution */
	virtual SCIP_DECL_CONSSEPASOL(scip_sepasol);

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
	virtual SCIP_RETCODE generate_Benders(
		SCIP *scip,
		SCIP_CONSHDLR *conshdlr,
		SCIP_SOL *sol,
		OsiCuts *cs);

	virtual SCIP_RETCODE sepaBenders(
		SCIP *scip,
		SCIP_CONSHDLR *conshdlr,
		SCIP_SOL *sol,
		SCIP_RESULT *result);

	virtual SCIP_RETCODE checkBenders(
		SCIP *scip,
		SCIP_CONSHDLR *conshdlr,
		SCIP_SOL *sol,
		SCIP_RESULT *result);

	/** generate Benders cuts */
	virtual void generateCuts(
		int size,  /**< [in] size of x */
		double *x, /**< [in] master solution */
		OsiCuts *cuts /**< [out] cuts generated */);

	/** generate Benders cuts */
	virtual void aggregateCuts(
			double ** cutvec, /**< [in] cut vector */
			double *  cutrhs, /**< [in] cut right-hand side */
			OsiCuts * cuts    /**< [out] cuts generated */);

	virtual void write_statistics();

protected:

	DecModel *  model_;            /**< DecModel object */
	BdSub *     bdsub_;            /**< pointer to cut generator */
	int         nvars_;            /**< number of original variables */
	SCIP_Var ** vars_;             /**< pointer array to original variables */
	int         naux_;             /**< number of auxiliary variables */
	double*     probability_;      /**< array of probability */

	/** simple statistics */
	vector<string> names_statistics_;
	unordered_map<string, int> count_statistics_;
	unordered_map<string, double> time_statistics_;
};

/** creates and captures a Benders constraint */
SCIP_RETCODE SCIPcreateConsBenders(
		SCIP * scip,
		SCIP_CONS ** cons,
		const char * name);

#endif /* SCIPCONSHDLRBENDERS_H_ */
