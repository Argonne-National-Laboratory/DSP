/*
 * DwWorker.h
 *
 *  Created on: Aug 27, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DANTZIGWOLFE_DWWORKER_H_
#define SRC_SOLVER_DANTZIGWOLFE_DWWORKER_H_

/** Coin */
#include "OsiSolverInterface.hpp"
/** Dsp */
#include "Model/DecModel.h"
#include "Utility/DspParams.h"
#include "Utility/DspMessage.h"

/**
 * This creates pricing subproblems. For each k, the subproblem is given by
 *   minimize  (c^k - pi^T H^k) x^k           - mu^k
 *   subject to
 *     lb^k <=              A^k x^k + B^k y^k <= ub^k,
 * where (pi,mu) are dual variables corresponding to the master rows.
 */
class DwWorker {
public:

	/** default constructor */
	DwWorker(DecModel * model, DspParams * par, DspMessage * message);

	/** default destructor */
	virtual ~DwWorker();

	/** generate variables */
	virtual DSP_RTN_CODE generateCols(
			int phase,                           /**< [in] phase of the master */
			const double* piA,                   /**< [in] piA */
			std::vector<int>& indices,           /**< [out] subproblem indices */
			std::vector<int>& statuses,          /**< [out] solution status */
			std::vector<double>& cxs,            /**< [out] solution times original objective coefficients */
			std::vector<double>& objs,           /**< [out] subproblem objective values */
			std::vector<CoinPackedVector*>& sols /**< [out] subproblem coupling column solutions */);

	/** get number of total subproblems */
	virtual int getNumSubprobs() {return nsubprobs_;}

	/** set column bounds */
	virtual void setColBounds(int size, const int* indices, const double* lbs, const double* ubs);

	/** set time limit */
	virtual void setTimeLimit(double limit);

	/** set gap tolerance */
	virtual void setGapTolerance(double gaptol);

private:

	/** set column bounds */
	virtual void setColBounds(int j, double lb, double ub);

protected:

	/** create subproblems */
	virtual DSP_RTN_CODE createSubproblems();

	/**
	 * This calculates and sets the objective coefficients for the subproblems.
	 */
	virtual DSP_RTN_CODE adjustObjFunction(
			int phase,        /**< [in] phase of the master */
			const double* piA /**< [in] dual variable times the constraint matrix */);

	/** solve subproblems */
	virtual DSP_RTN_CODE solveSubproblems();

	/** reset subproblems */
	virtual DSP_RTN_CODE resetSubproblems();

public:

	DecModel * model_;     /**< DecModel object */
	DspParams * par_;      /**< parameters */
	DspMessage * message_; /**< message */

protected:

	OsiSolverInterface** si_; /**< solver interface */

	//DwSub* sub_; /**< subproblem solver */
	double** sub_objs_; /**< subproblem objective coefficients */
	double** sub_clbd_; /**< subproblem column lower bounds */
	double** sub_cubd_; /**< subproblem column lower bounds */
	bool** coupled_; /**< indicate coupling columns (size of parProcIdxSize_ of getNumSubproblemCouplingCols(s)) */

	int parProcIdxSize_; /**< number of subproblems for this worker */
	int * parProcIdx_;   /**< subproblem indices for this worker */

	int nsubprobs_; /**< number of total subproblems */
};

#endif /* SRC_SOLVER_DANTZIGWOLFE_DWWORKER_H_ */
