/*
 * SCIPconshdlrIntBenders.h
 */

#ifndef SCIPCONSHDLRINTBENDERS_H_
#define SCIPCONSHDLRINTBENDERS_H_

#include "Solver/Benders/SCIPconshdlrDrBenders.h"

/** A base class for implementing integer Benders constraint handler */
class SCIPconshdlrIntBenders : public SCIPconshdlrDrBenders
{
public:
	/** default constructor */
	SCIPconshdlrIntBenders(SCIP *scip, const char *name, int sepapriority, int param_seps_solver);

	/** default constructor */
	virtual ~SCIPconshdlrIntBenders();

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

protected:
	virtual double compute_weighted_sum(
		const double *recourse_values /**< [in] recourse values */
	);

	virtual double compute_approximate_recourse(
		SCIP *scip,	   /**< [in] scip pointer */
		SCIP_SOL *sol, /**< [in] solution to evaluate */
		const double *recourse_values /**< [in] recourse values */
	);
	
	virtual SCIP_RETCODE tryPrimalSolution(
		SCIP *scip,	   /**< [in] scip pointer */
		SCIP_SOL *sol, /**< [in] solution to evaluate */
		const double *recourse_values /**< [in] recourse values */);

	virtual SCIP_RETCODE addNoGoodCut(
		SCIP *scip,				 /**< [in] scip pointer */
		SCIP_CONSHDLR *conshdlr, /**< [in] constraint handler that creates the row */
		SCIP_RESULT *result /**< [out] result */);

	virtual SCIP_RETCODE tryIntOptimalityCut(
		SCIP *scip,					   /**< [in] scip pointer */
		SCIP_CONSHDLR *conshdlr,	   /**< [in] constraint handler that creates the row */
		SCIP_SOL *sol,				   /**< [in] solution to evaluate */
		const double *recourse_values, /**< [in] recourse values */
		SCIP_RESULT *result /**< [out] result */);

	virtual SCIP_RETCODE addIntOptimalityCut(
		SCIP *scip,				 /**< [in] scip pointer */
		SCIP_CONSHDLR *conshdlr, /**< [in] constraint handler that creates the row */
		double exact_recourse,	 /**< [in] exact recourse value */
		SCIP_RESULT *result /**< [out] result */);

	virtual bool check_binary_solution_pool(
		SCIP *scip, /**< [in] scip pointer */
		SCIP_SOL *sol /**< [in] solution to evaluate */);

	virtual void add_binary_solution_pool(
		SCIP *scip, /**< [in] scip pointer */
		SCIP_SOL *sol /**< [in] solution to evaluate */);

protected:
	vector<vector<int>> binary_solution_pool_; /**< binary solution pool */
};

/** creates and captures a Benders constraint */
SCIP_RETCODE SCIPcreateConsIntBenders(
	SCIP *scip,
	SCIP_CONS **cons,
	const char *name);

#endif /* SCIPCONSHDLRINTBENDERS_H_ */
