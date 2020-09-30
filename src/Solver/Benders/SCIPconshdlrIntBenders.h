/*
 * SCIPconshdlrIntBenders.h
 */

#ifndef SCIPCONSHDLRINTBENDERS_H_
#define SCIPCONSHDLRINTBENDERS_H_

#include "Solver/Benders/SCIPconshdlrBenders.h"

/** A base class for implementing Benders constraint handler */
class SCIPconshdlrIntBenders : public SCIPconshdlrBenders
{
public:
	/** default constructor */
	SCIPconshdlrIntBenders(SCIP *scip, const char *name, int sepapriority);

	/** default constructor */
	virtual ~SCIPconshdlrIntBenders();

	/** constraint enforcing method of constraint handler for LP solutions */
	virtual SCIP_DECL_CONSENFOLP(scip_enfolp);

	/** constraint enforcing method of constraint handler for pseudo solutions */
	virtual SCIP_DECL_CONSENFOPS(scip_enfops);

	/** feasibility check method of constraint handler for primal solutions */
	virtual SCIP_DECL_CONSCHECK(scip_check);

protected:
	/** evaluate recourse: 
	 * This function is called only when the recourse has integer variables.
	*/
	virtual SCIP_RETCODE evaluateRecourse(
		SCIP *scip,	   /**< [in] scip pointer */
		SCIP_SOL *sol, /**< [in] solution to evaluate */
		double *values /**< [out] evaluated recourse values */);

	virtual SCIP_RETCODE addNoGoodCut(
		SCIP *scip,				 /**< [in] scip pointer */
		SCIP_CONSHDLR *conshdlr, /**< [in] constraint handler that creates the row */
		SCIP_RESULT *result /**< [out] result */);

	virtual SCIP_RETCODE addIntOptimalityCut(
		SCIP *scip,				 /**< [in] scip pointer */
		SCIP_CONSHDLR *conshdlr, /**< [in] constraint handler that creates the row */
		double exact_recourse,	 /**< [in] exact recourse value */
		double recourse_lb,		 /**< [in] lower bound of recourse value */
		SCIP_RESULT *result /**< [out] result */);

	virtual bool check_binary_solution_pool(
		SCIP *scip,	   /**< [in] scip pointer */
		SCIP_SOL *sol, /**< [in] solution to evaluate */
		bool append = false /**< [in] whether sol is appended to the pool */);

protected:
	vector<vector<int>> binary_solution_pool_; /**< binary solution pool */
};

/** creates and captures a Benders constraint */
SCIP_RETCODE SCIPcreateConsIntBenders(
	SCIP *scip,
	SCIP_CONS **cons,
	const char *name);

#endif /* SCIPCONSHDLRINTBENDERS_H_ */
