/*
 * StoParam.h
 *
 *  Created on: Nov 6, 2014
 *      Author: kibaekkim
 */

#ifndef STOPARAM_H_
#define STOPARAM_H_

#include <limits>
#include "Utility/StoMacros.h"

enum TssDdMasterSolver
{
	Simplex = 0,
	IPM,
	IPM_Feasible,
	DSBM, /**< doubly stabilized bundle method */
	Subgradient
};

class StoParam
{
public:
	StoParam():
		logLevel_(0),
		numCores_(1),
		TssBdNumAugScenarios_(0),
		TssBdAugScenarios_(NULL),
		TssBdBendersPriority_(-2000000),
		TssDdNumProcIdx_(0),
		TssDdProcIdxSet_(NULL),
		TssDdMasterSolver_(IPM_Feasible),
		TssDdMasterNumCutsPerIter_(1),
		TssDdEnableTrustRegion_(1),
		TssDdTrustRegionSize_(100.),
		TssDdDisableTrustRegionDecrease_(false),
		TssDdCacheRecourse_(0),
		TssDdAddFeasCuts_(-1),
		TssDdAddOptCuts_(-1),
		TssDdEvalUb_(1),
		TssDdDualVarsLog_(0),
		TssDdStoppingTol_(1.e-5),
		ScipDisplayFreq_(100),
		ScipLimitsGap_(0.0),
		ScipLimitsTime_(1.0e+20)
	{
		nodeLimit_ = -1;
		iterLimit_ = std::numeric_limits<int>::max();
		wtimeLimit_ = std::numeric_limits<double>::max();
		relaxIntegrality_[0] = false;
		relaxIntegrality_[1] = false;
		relaxIntegralityAll_ = false;
	}
	virtual ~StoParam()
	{
		FREE_ARRAY_PTR(TssBdAugScenarios_);
		FREE_ARRAY_PTR(TssDdProcIdxSet_);
	}

public:

	/** print level:
	 * 0: nothing
	 * 1: minimum
	 *    Bd, Dd - iteration information,
	 *             nothing from optimization solvers (e.g., OOQP, SCIP, Clp)
	 * De - set as solver print level
	 * */
	int logLevel_;

	/** number of cores for OpenMP solver (Bd only) */
	int numCores_;

	/** Branch-and-cut node limit:
	 * De - solver node limit
	 * Bd - solver node limit for Phase 2
	 * TODO Dd - solver node limit for subproblems and feasibility recovery
	 * */
	int nodeLimit_;

	/** iteration limit:
	 * De - solver limit
	 * Bd - L-shaped iteration limit for LP
	 * Dd - dual decomposition iteration limit
	 */
	int iterLimit_;

	/** wall clock limit */
	double wtimeLimit_;

	bool relaxIntegrality_[2];

	bool relaxIntegralityAll_; /* TODO: Not implemented in TssBd */

	/**
	 * Benders decomposition
	 * */

	/** Benders augmented scenarios */
	int TssBdNumAugScenarios_;
	int * TssBdAugScenarios_;

	/** Benders cut priority (refer CONSHDLR_SEPAPRIORITY of SCIP constraint handler) */
	int TssBdBendersPriority_;


	/**
	 * Dual decomposition
	 * */

	/** set of scenarios for the current process */
	int TssDdNumProcIdx_;
	int * TssDdProcIdxSet_;

	/** Lagrangian master problem solver (see TssDdMasterSolver)*/
	int TssDdMasterSolver_;

	/** number of cuts per iteration */
	int TssDdMasterNumCutsPerIter_;

	/** Dual decomposition enables trust region if 1 (default); 0 otherwise */
	int TssDdEnableTrustRegion_;

	/** Initial trust region size */
	double TssDdTrustRegionSize_;

	/** Disable decrease in trust region size */
	bool TssDdDisableTrustRegionDecrease_;

	/** cache recourse problems to get upper bounds if 1; 1 otherwise (default) */
	int TssDdCacheRecourse_;

	/** add feasibility cuts:
	 * < 0: do not add feasibility cut
	 *   0: do only at the first iteration
	 * > 0: do every TssDdAddFeasCuts_ iterations
	 */
	int TssDdAddFeasCuts_;

	/** add optimality cuts:
	 * < 0: do not add optimality cut
	 *   0: do only at the first iteration
	 * > 0: do every TssDdAddOptCuts_ iterations
	 */
	int TssDdAddOptCuts_;

	/** evaluate upper bound:
	 * < 0: do not evaluate upper bound
	 *   0: do only at the first iteration
	 * > 0: do every TssDdEvalUb_ iterations
	 */
	int TssDdEvalUb_;

	/** Store changes of distance of dual variables (only for experiments) */
	int TssDdDualVarsLog_;

	/** Stopping gap tolerance
	 *   - relative gap between master objective and dual bound
	 *   - relative gap between primal bound and dual bound
	 */
	double TssDdStoppingTol_;

	int ScipDisplayFreq_;
	double ScipLimitsGap_;

	/** time limit for SCIP solver */
	double ScipLimitsTime_;
};

#endif /* STOPARAM_H_ */
