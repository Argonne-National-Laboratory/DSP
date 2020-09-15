#define DSP_DEBUG

#include "Utility/DspMessage.h"
#include "Model/TssModel.h"
#include "Solver/Benders/SCIPheurIntBenders.h"

SCIP_DECL_HEUREXEC(SCIPheurIntBenders::scip_exec) {
    *result = SCIP_DIDNOTRUN;
    DSPdebugMessage("node infeasible? %s\n", nodeinfeasible ? "yes" : "no");

    /* only call heuristic, if an optimal LP solution is at hand or if relaxation solution is available */
    if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL && ! SCIPisRelaxSolValid(scip) )
        return SCIP_OKAY;
 
    /* only call heuristic, if the LP objective value is smaller than the cutoff bound */
    if( SCIPisGE(scip, SCIPgetLPObjval(scip), SCIPgetCutoffbound(scip)) )
        return SCIP_OKAY;

    int nvars = SCIPgetNOrigVars(scip);
    SCIP_VAR** vars = SCIPgetOrigVars(scip);
    SCIP_Real* vals = NULL; /**< variable values */

    /** allocate memory */
    SCIP_CALL(SCIPallocMemoryArray(scip, &vals, nvars));

	/** get current solution */
	SCIP_CALL(SCIPgetSolVals(scip, NULL, nvars, vars, vals));

    /** is the solution feasible integral? */
    bool is_feas_integral = true;
    for (int j = 0; j < nvars; ++j) {
        if (SCIPvarGetType(vars[j]) != SCIP_VARTYPE_CONTINUOUS) {
            if (!SCIPisFeasIntegral(scip, vals[j])) {
                is_feas_integral = false;
                break;
            }
        }
    }
    DSPdebugMessage("is feasible integral? %s\n", is_feas_integral ? "yes" : "no");

    /* only call heuristic, if an optimal LP solution is integer feasible */
    if (is_feas_integral) {
        // TODO:
    }

	/** free memory */
	SCIPfreeMemoryArray(scip, &vals);

    return SCIP_OKAY;
}