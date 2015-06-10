/*
 * StoRtnCodes.h
 *
 *  Created on: Sep 25, 2014
 *      Author: kibaekkim
 */

#ifndef STORTNCODES_H_
#define STORTNCODES_H_

/*
 * This defines return codes.
 */

typedef int STO_RTN_CODE;

#define STO_RTN_OK    0
#define STO_RTN_ERR   1

#define STO_STAT_OPTIMAL            3000
#define STO_STAT_PRIM_INFEASIBLE    3001
#define STO_STAT_DUAL_INFEASIBLE    3002
#define STO_STAT_LIM_ITERorTIME     3004
#define STO_STAT_STOPPED_GAP        3005
#define STO_STAT_STOPPED_NODE       3006
#define STO_STAT_STOPPED_TIME       3007
#define STO_STAT_STOPPED_USER       3008
#define STO_STAT_STOPPED_SOLUTION   3009
#define STO_STAT_STOPPED_ITER       3010
#define STO_STAT_STOPPED_UNKNOWN    3011
#define STO_STAT_STOPPED_MPI        3012
#define STO_STAT_ABORT              3013
#define STO_STAT_LIM_PRIM_OBJ       3014
#define STO_STAT_LIM_DUAL_OBJ       3015
#define STO_STAT_FEASIBLE           3016
#define STO_STAT_UNKNOWN            3999

#define STO_RTN_MSG_BODY "Error code %d in %s::%s"

#define STO_RTN_CHECK(__Rtn,__Func,__Class)                       \
	{                                                             \
		STO_RTN_CODE __rtn = __Rtn;                               \
		if (__rtn != STO_RTN_OK) {                                \
			printf(STO_RTN_MSG_BODY"\n", __rtn, __Func, __Class); \
		}                                                         \
	}

#define STO_RTN_CHECK_RTN(__Rtn,__Func,__Class)                   \
	{                                                             \
		STO_RTN_CODE __rtn = __Rtn;                               \
		if (__rtn != STO_RTN_OK) {                                \
			printf(STO_RTN_MSG_BODY"\n", __rtn, __Func, __Class); \
			return;                                               \
		}                                                         \
	}

#define STO_RTN_CHECK_RTN_CODE(__Rtn,__Func,__Class)              \
	{                                                             \
		STO_RTN_CODE __rtn = __Rtn;                               \
		if (__rtn != STO_RTN_OK) {                                \
			printf(STO_RTN_MSG_BODY"\n", __rtn, __Func, __Class); \
			return __rtn;                                         \
		}                                                         \
	}

#define STO_RTN_CHECK_THROW(__Rtn,__Func,__Class)                        \
	{                                                                    \
		STO_RTN_CODE __rtn = __Rtn;                                      \
		if (__rtn != STO_RTN_OK) {                                       \
			char __tmpstr[128];                                          \
			sprintf(__tmpstr, STO_RTN_MSG_BODY, __rtn, __Func, __Class); \
			throw __tmpstr;                                              \
		}                                                                \
	}

#endif /* STORTNCODES_H_ */
