/**
 * DspOsiScip.h
 *
 * 12/12/2019
 * Kibaek Kim
 */
 
#ifndef SRC_SOLVERINTERFACE_DSPOSISCIP_H_
#define SRC_SOLVERINTERFACE_DSPOSISCIP_H_

#ifdef DSP_HAS_SCIP

#include "SolverInterface/OsiScipSolverInterface.hpp"
#include "SolverInterface/SCIPbranchruleLB.h"

inline void convertOsiScipToDspStatus(OsiScipSolverInterface* si, int& status) {
    int scipstat = SCIPgetStatus(si->getScip());
    switch(scipstat) {
    case SCIP_STATUS_USERINTERRUPT:
        status = DSP_STAT_STOPPED_USER;
        break;
    case SCIP_STATUS_NODELIMIT:
    case SCIP_STATUS_TOTALNODELIMIT:
        status = DSP_STAT_STOPPED_NODE;
        break;
    case SCIP_STATUS_TIMELIMIT:
        status = DSP_STAT_STOPPED_TIME;
        break;
    case SCIP_STATUS_GAPLIMIT:
        status = DSP_STAT_STOPPED_GAP;
        break;
    case SCIP_STATUS_SOLLIMIT:
    case SCIP_STATUS_BESTSOLLIMIT:
        status = DSP_STAT_STOPPED_SOLUTION;
        break;
    case SCIP_STATUS_STALLNODELIMIT:
    case SCIP_STATUS_MEMLIMIT:
    case SCIP_STATUS_RESTARTLIMIT:
        status = DSP_STAT_STOPPED_UNKNOWN;
        break;
    case SCIP_STATUS_OPTIMAL:
        status = DSP_STAT_OPTIMAL;
        break;
    case SCIP_STATUS_INFEASIBLE:
        status = DSP_STAT_PRIM_INFEASIBLE;
        break;
    case SCIP_STATUS_UNBOUNDED:
    case SCIP_STATUS_INFORUNBD:
        status = DSP_STAT_DUAL_INFEASIBLE;
        break;
    case SCIP_STATUS_UNKNOWN:
    default:
        status = DSP_STAT_UNKNOWN;
        break;
    }
}

#endif /* DSP_HAS_SCIP */

#endif /* SRC_SOLVERINTERFACE_DSPOSISCIP_H_ */