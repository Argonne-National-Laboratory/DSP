/**
 * DspOsiCpx.h
 *
 * 12/12/2019
 * Kibaek Kim
 */
 
#ifndef SRC_SOLVERINTERFACE_DSPOSICPX_H_
#define SRC_SOLVERINTERFACE_DSPOSICPX_H_

#ifdef DSP_HAS_CPX

#include "cplex.h"
#include "OsiCpxSolverInterface.hpp"

inline void convertOsiCpxToDspStatus(OsiCpxSolverInterface* si, int& status) {
	int probtype = CPXgetprobtype(si->getEnvironmentPtr(), si->getLpPtr());
	int stat = CPXgetstat(si->getEnvironmentPtr(), si->getLpPtr());
	if (probtype == CPXPROB_LP) {
		switch(stat) {
		case CPX_STAT_OPTIMAL:
		case CPX_STAT_OPTIMAL_INFEAS:
			status = DSP_STAT_OPTIMAL;
			break;
		case CPX_STAT_INFEASIBLE:
			status = DSP_STAT_PRIM_INFEASIBLE;
			break;
		case CPX_STAT_UNBOUNDED:
			status = DSP_STAT_DUAL_INFEASIBLE;
			break;
		case CPX_STAT_ABORT_OBJ_LIM:
		case CPX_STAT_ABORT_PRIM_OBJ_LIM:
			status = DSP_STAT_LIM_PRIM_OBJ;
			break;
		case CPX_STAT_ABORT_DUAL_OBJ_LIM:
			status = DSP_STAT_LIM_DUAL_OBJ;
			break;
		case CPX_STAT_ABORT_IT_LIM:
			status = DSP_STAT_STOPPED_ITER;
			break;
		case CPX_STAT_ABORT_TIME_LIM:
			status = DSP_STAT_STOPPED_TIME;
			break;
		case CPX_STAT_NUM_BEST:
		case CPX_STAT_FEASIBLE:
			status = DSP_STAT_FEASIBLE;
			break;
		case CPX_STAT_ABORT_USER:
			status = DSP_STAT_STOPPED_USER;
			break;
		default:
			status = DSP_STAT_UNKNOWN;
			break;
		}
	} else if (probtype == CPXPROB_MILP) {
		switch(stat) {
		case CPXMIP_OPTIMAL:
		case CPXMIP_OPTIMAL_TOL:
		case CPXMIP_OPTIMAL_INFEAS:
			status = DSP_STAT_OPTIMAL;
			break;
		case CPXMIP_INFEASIBLE:
		case CPXMIP_NODE_LIM_INFEAS:
		case CPXMIP_TIME_LIM_INFEAS:
			status = DSP_STAT_PRIM_INFEASIBLE;
			break;
		case CPXMIP_UNBOUNDED:
		case CPXMIP_INForUNBD:
			status = DSP_STAT_DUAL_INFEASIBLE;
			break;
		case CPXMIP_NODE_LIM_FEAS:
			status = DSP_STAT_STOPPED_NODE;
			break;
		case CPXMIP_TIME_LIM_FEAS:
			status = DSP_STAT_STOPPED_TIME;
			break;
		case CPXMIP_FAIL_FEAS:
		case CPXMIP_FAIL_INFEAS:
		case CPXMIP_MEM_LIM_FEAS:
		case CPXMIP_MEM_LIM_INFEAS:
		case CPXMIP_ABORT_FEAS:
		case CPXMIP_ABORT_INFEAS:
		case CPXMIP_FAIL_FEAS_NO_TREE:
		case CPXMIP_FAIL_INFEAS_NO_TREE:
			status = DSP_STAT_ABORT;
			break;
		default:
			status = DSP_STAT_UNKNOWN;
			break;
		}
	} else {
		status = DSP_STAT_UNKNOWN;
	}
}

#endif

#endif /* SRC_SOLVERINTERFACE_DSPOSICPX_H_ */