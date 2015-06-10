/*
 * SCIPbranchruleLB.cpp
 *
 *  Created on: Apr 9, 2015
 *      Author: kibaekkim
 */

#include <Solver/SCIPbranchruleLB.h>


/** branching execution method for fractional LP solutions
 *
 *  @see SCIP_DECL_BRANCHEXECLP(x) in @ref type_branch.h
 */
SCIP_DECL_BRANCHEXECLP(SCIPbranchruleLB::scip_execlp)
{
	checkLowerBound(scip, result);
	return SCIP_OKAY;
}

/** branching execution method for external candidates
 *
 *  @see SCIP_DECL_BRANCHEXECEXT(x) in @ref type_branch.h
 */
SCIP_DECL_BRANCHEXECEXT(SCIPbranchruleLB::scip_execext)
{
	checkLowerBound(scip, result);
	return SCIP_OKAY;
}

/** branching execution method for not completely fixed pseudo solutions
 *
 *  @see SCIP_DECL_BRANCHEXECPS(x) in @ref type_branch.h
 */
SCIP_DECL_BRANCHEXECPS(SCIPbranchruleLB::scip_execps)
{
	checkLowerBound(scip, result);
	return SCIP_OKAY;
}

/** check lower bound */
void SCIPbranchruleLB::checkLowerBound(SCIP * scip, SCIP_RESULT * result)
{
	*result = SCIP_DIDNOTRUN;

	if (lowerbound_ <= -COIN_DBL_MAX) return;

	double objval = SCIPgetLPObjval(scip);

	switch(SCIPgetLPSolstat(scip))
	{
	case SCIP_LPSOLSTAT_OPTIMAL:
	case SCIP_LPSOLSTAT_OBJLIMIT:
	case SCIP_LPSOLSTAT_ITERLIMIT:
	case SCIP_LPSOLSTAT_TIMELIMIT:
		if (SCIPisLT(scip, lowerbound_, objval))
		{
			printf("-> branchrule: lowerbound_ %e objval %e\n", lowerbound_, objval);
			*result = SCIP_CUTOFF;
		}
		break;
	case SCIP_LPSOLSTAT_NOTSOLVED:
	case SCIP_LPSOLSTAT_INFEASIBLE:
	case SCIP_LPSOLSTAT_UNBOUNDEDRAY:
	case SCIP_LPSOLSTAT_ERROR:
		break;
	}
}
