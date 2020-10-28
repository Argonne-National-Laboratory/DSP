/*
 * SCIPconshdlrDrBenders.h
 */

#ifndef SCIPCONSHDLRDRBENDERS_H_
#define SCIPCONSHDLRDRBENDERS_H_

#include "Solver/Benders/SCIPconshdlrBenders.h"

/** A base class for implementing distributionally robust Benders constraint handler */
class SCIPconshdlrDrBenders : public SCIPconshdlrBenders
{
public:
	/** default constructor */
	SCIPconshdlrDrBenders(SCIP *scip, const char *name, int sepapriority, int param_seps_solver);

	/** default constructor */
	virtual ~SCIPconshdlrDrBenders();

	/** constraint enforcing method of constraint handler for LP solutions */
	virtual SCIP_DECL_CONSENFOLP(scip_enfolp);

	/** constraint enforcing method of constraint handler for pseudo solutions */
	virtual SCIP_DECL_CONSENFOPS(scip_enfops);

	/** feasibility check method of constraint handler for primal solutions */
	virtual SCIP_DECL_CONSCHECK(scip_check);

	/** separation method of constraint handler for LP solution */
	virtual SCIP_DECL_CONSSEPALP(scip_sepalp);

	/** separation method of constraint handler for arbitrary primal solution */
	virtual SCIP_DECL_CONSSEPASOL(scip_sepasol);

	/** set model pointer */
	virtual void setDecModel(DecModel *model);

protected:
	/** evaluate recourse functions */
	virtual SCIP_RETCODE evaluateRecourse(
		SCIP *scip,	   /**< [in] scip pointer */
		SCIP_SOL *sol, /**< [in] solution to evaluate */
		double *values /**< [out] evaluated recourse values */);

	/** create distribution separation problem */
	virtual void create_distsepa_problem();

	/** compuate probability */
	virtual SCIP_RETCODE computeProbability(
		const double *recourse /**< [in] recourse values */);

	/** distributionally robust Benders separation */
	virtual SCIP_RETCODE sepaDrBenders(
		SCIP *scip, /**< [in] scip pointer */
		SCIP_CONSHDLR *conshdlr,
		SCIP_SOL *sol, /**< [in] solution to evaluate */
		SCIP_RESULT *result /**< [out] result */);

	/** distributionally robust Benders separation (checking only) */
	virtual SCIP_RETCODE checkDrBenders(
		SCIP *scip, /**< [in] scip pointer */
		SCIP_CONSHDLR *conshdlr,
		SCIP_SOL *sol, /**< [in] solution to evaluate */
		SCIP_RESULT *result /**< [out] result */);

	DspOsi *drosi_; /**< distribution separation problem */
};

/** creates and captures a Benders constraint */
SCIP_RETCODE SCIPcreateConsDrBenders(
	SCIP *scip,
	SCIP_CONS **cons,
	const char *name);

#endif /* SCIPCONSHDLRDRBENDERS_H_ */
