#ifndef SRC_SOLVERINTERFACE_DSPOSIGRB_H_
#define SRC_SOLVERINTERFACE_DSPOSIGRB_H_

#ifdef DSP_HAS_GRB

#include "gurobi_c.h"
#include "OsiGrbSolverInterface.hpp"

inline void convertOsiGrbToDspStatus(OsiGrbSolverInterface* si, int& status){
    int Grbstatus;
    GRBgetintattr(si->getLpPtr(), GRB_INT_ATTR_STATUS, &Grbstatus);
    switch(Grbstatus){
        case GRB_LOADED:
        case GRB_OPTIMAL:
            status = DSP_STAT_OPTIMAL;
            break;
        case GRB_INFEASIBLE:
            status = DSP_STAT_PRIM_INFEASIBLE;
            break;
        //case GRB_INF_OR_UNBD:
        case GRB_UNBOUNDED:
            status = DSP_STAT_DUAL_INFEASIBLE;
            break;
        case GRB_CUTOFF:
        case GRB_ITERATION_LIMIT:
            status = DSP_STAT_STOPPED_ITER;
			break;
        case GRB_NODE_LIMIT:
            status = DSP_STAT_STOPPED_NODE;
            break;
        case GRB_TIME_LIMIT:
            status = DSP_STAT_STOPPED_TIME;
            break;
        case GRB_SOLUTION_LIMIT:
        case GRB_INTERRUPTED:
            status = DSP_STAT_STOPPED_USER;
            break;
        case GRB_NUMERIC:
        case GRB_SUBOPTIMAL:
            status = DSP_STAT_STOPPED_SOLUTION;
            break;
        case GRB_INPROGRESS:
        case GRB_USER_OBJ_LIMIT:
        default:
            status = DSP_STAT_UNKNOWN;
            break;
    }
}

#endif

#endif //SRC_SOLVERINTERFACE_DSPOSIGRB_H_
