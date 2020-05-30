/*
 * DspRtnCodes.h
 *
 *  Refactored on: April 8, 2016
 *  Created on: Sep 25, 2014
 *      Author: kibaekkim
 */

#ifndef DSPRTNCODES_H_
#define DSPRTNCODES_H_

/*
 * This defines return codes.
 */

typedef int DSP_RTN_CODE;

#define DSP_RTN_OK    0
#define DSP_RTN_ERR   1

#define DSP_STAT_OPTIMAL            3000
#define DSP_STAT_PRIM_INFEASIBLE    3001
#define DSP_STAT_DUAL_INFEASIBLE    3002
#define DSP_STAT_LIM_ITERorTIME     3004
#define DSP_STAT_STOPPED_GAP        3005
#define DSP_STAT_STOPPED_NODE       3006
#define DSP_STAT_STOPPED_TIME       3007
#define DSP_STAT_STOPPED_USER       3008
#define DSP_STAT_STOPPED_SOLUTION   3009
#define DSP_STAT_STOPPED_ITER       3010
#define DSP_STAT_STOPPED_UNKNOWN    3011
#define DSP_STAT_STOPPED_MPI        3012
#define DSP_STAT_ABORT              3013
#define DSP_STAT_LIM_PRIM_OBJ       3014
#define DSP_STAT_LIM_DUAL_OBJ       3015
#define DSP_STAT_FEASIBLE           3016
#define DSP_STAT_RUN_HEURISTICS     3017
#define DSP_STAT_LIM_INFEAS         3018
#define DSP_STAT_STOPPED_NUMERICS   3019
#define DSP_STAT_MW_STOP            3100 /**< stop signal for master-worker framework */
#define DSP_STAT_MW_CONTINUE        3101 /**< continue signal for master-worker framework */
#define DSP_STAT_MW_EXACT           3102 /**< force to evaluate exactly */
#define DSP_STAT_MW_RESOLVE         3103 /**< need to resolve */
#define DSP_STAT_NOT_SOLVED         3998
#define DSP_STAT_UNKNOWN            3999

#define DSP_RTN_MSG_BODY "Error code %d in %s:%d"

#define DSP_RTN_CHECK(__Rtn)                       \
	{                                                             \
		DSP_RTN_CODE __rtn = __Rtn;                               \
		if (__rtn != DSP_RTN_OK) {                                \
			printf(DSP_RTN_MSG_BODY"\n", __rtn, __FILE__, __LINE__); \
		}                                                         \
	}

#define DSP_RTN_CHECK_RTN(__Rtn)                   \
	{                                                             \
		DSP_RTN_CODE __rtn = __Rtn;                               \
		if (__rtn != DSP_RTN_OK) {                                \
			printf(DSP_RTN_MSG_BODY"\n", __rtn, __FILE__, __LINE__); \
			return;                                               \
		}                                                         \
	}

#define DSP_RTN_CHECK_RTN_CODE(__Rtn)              \
	{                                                             \
		DSP_RTN_CODE __rtn = __Rtn;                               \
		if (__rtn != DSP_RTN_OK) {                                \
			printf(DSP_RTN_MSG_BODY"\n", __rtn, __FILE__, __LINE__); \
			return __rtn;                                         \
		}                                                         \
	}

#define DSP_RTN_CHECK_THROW(__Rtn)                        \
	{                                                                    \
		DSP_RTN_CODE __rtn = __Rtn;                                      \
		if (__rtn != DSP_RTN_OK) {                                       \
			char __tmpstr[128];                                          \
			sprintf(__tmpstr, DSP_RTN_MSG_BODY, __rtn, __FILE__, __LINE__); \
			throw __tmpstr;                                              \
		}                                                                \
	}

#endif /* DSPRTNCODES_H_ */
