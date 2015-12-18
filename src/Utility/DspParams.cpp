/*
 * DspParams.cpp
 *
 *  Created on: Oct 20, 2015
 *      Author: kibaekkim
 */

#include <Utility/DspParams.h>

#define MAX_INT_NUM numeric_limits<int>::max()
#define MAX_DBL_NUM numeric_limits<double>::max()

DspParams::DspParams()
{
	initBoolParams();
	initIntParams();
	initDblParams();
	initStrParams();
	initBoolPtrParams();
	initIntPtrParams();
}

DspParams::~DspParams()
{
	BoolPtrParams_.deleteParam("RELAX_INTEGRALITY");
	IntPtrParams_.deleteParam("BD/ARR_AUG_SCENS");
	IntPtrParams_.deleteParam("ARR_PROC_IDX");
}

void DspParams::initBoolParams()
{
	/** enable trust region */
	BoolParams_.createParam("DD/TR", true);

	/** enable decreasing trust region */
	BoolParams_.createParam("DD/TR/DECREASE", true);

	/** cache recourse models */
	BoolParams_.createParam("DD/CACHE_RECOURSE", true);

	/** log dual variable values */
	BoolParams_.createParam("DD/LOG_DUAL_VARS", false);
}

void DspParams::initIntParams()
{
	/** print level */
	IntParams_.createParam("LOG_LEVEL", 1);

	/** branch-and-cut node limit */
	IntParams_.createParam("NODE_LIM", MAX_INT_NUM);

	/** iteration limit */
	IntParams_.createParam("ITER_LIM", MAX_INT_NUM);

	/** number of cores used in OpenMP library (Benders only) */
	IntParams_.createParam("BD/NUM_CORES", 1);

	/** Benders cut priority (refer CONSHDLR_SEPAPRIORITY of SCIP constraint handler */
	IntParams_.createParam("BD/CUT_PRIORITY", -200000);

	/** Benders lower bound methods:
	 * 0 = solve separate LP relaxation problems;
	 * 1 = solve separate MILP relaxation problems */
	IntParams_.createParam("BD/INIT_LB_ALGO", SEPARATE_LP);

	/** algorithm for the master */
	IntParams_.createParam("DD/MASTER_ALGO", IPM_Feasible);

	/** number of cuts to the master per iteration */
	IntParams_.createParam("DD/NUM_CUTS_PER_ITER", 1);

	/** add feasibility cuts */
	IntParams_.createParam("DD/FEAS_CUTS", -1);

	/** add optimality cuts */
	IntParams_.createParam("DD/OPT_CUTS", -1);

	/** evaluate upper bound */
	IntParams_.createParam("DD/EVAL_UB", -1);

	/** display frequency */
	IntParams_.createParam("SCIP/DISPLAY_FREQ", 100);
}

void DspParams::initDblParams()
{
	/** wall clock limit */
	DblParams_.createParam("WALL_LIM", MAX_DBL_NUM);

	/** initial trust region size */
	DblParams_.createParam("DD/TR/SIZE", 100);

	/** stopping tolerance */
	DblParams_.createParam("DD/STOP_TOL", 1.0e-5);

	/** branch-and-bound gap tolerance */
	DblParams_.createParam("SCIP/GAP_TOL", 0.0);

	/** time limit */
	DblParams_.createParam("SCIP/TIME_LIM", MAX_DBL_NUM);
}

void DspParams::initStrParams()
{
	/** no string parameters */
}

void DspParams::initBoolPtrParams()
{
	/** array of indicating relaxation of integer variables;
	 * each element represents a stage. */
	BoolPtrParams_.createParam("RELAX_INTEGRALITY", 2);
	BoolPtrParams_.setParam("RELAX_INTEGRALITY", 0, false);
	BoolPtrParams_.setParam("RELAX_INTEGRALITY", 1, false);
}

void DspParams::initIntPtrParams()
{
	/** array of scenarios for the current process */
	IntPtrParams_.createParam("ARR_PROC_IDX");

	/** array of augmented scenarios */
	IntPtrParams_.createParam("BD/ARR_AUG_SCENS");
}
