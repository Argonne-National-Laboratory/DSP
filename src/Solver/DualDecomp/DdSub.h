/*
 * DdSub.h
 *
 *  Renamed on: Feb 16, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DDSUB_H_
#define SRC_SOLVER_DDSUB_H_

#include <typeinfo>

/** DSP */
#include "Model/DecModel.h"
#include "Solver/DecSolver.h"
#include "SolverInterface/SolverInterface.h"

/**
 * This defines a class for solving a dual decomposition subproblem.
 */
class DdSub : public DecSolver
{
public:

	/** default constructor */
	DdSub(int s, DspParams * par, DecModel * model, DspMessage * message);

	/** default destructor */
	virtual ~DdSub();

	/** initialize */
	virtual DSP_RTN_CODE init();

	/** finalize */
	virtual DSP_RTN_CODE finalize() {return DSP_RTN_OK;}

	/** solve */
	virtual DSP_RTN_CODE solve();

private:

	/** create problem */
	DSP_RTN_CODE createProblem();

	/** add cut generator */
	DSP_RTN_CODE addCutGenerator();

public:

	/** update problem */
	DSP_RTN_CODE updateProblem(
			double * lambda,
			double primal_bound = COIN_DBL_MAX);

	/** push cuts */
	DSP_RTN_CODE pushCuts(OsiCuts * cuts);

	/** set wall clock time limit */
	void setTimeLimit(double sec);

	/** set accuracy tolerance */
	void setGapTol(double tol);

	/** set print level */
	void setPrintLevel(int level);

public:

	SolverInterface * si_; /**< solver interface */

	int sind_;           /**< scenario index */
	int nrows_coupling_; /**< number of coupling constraints for the subproblem (dimension of lambda) */
	int ncols_coupling_; /**< number of coupling variables for the subproblem
						   *  (dimension of the part of the solution returned to the master problem) */

	double theta_;  /**< theta */
	double gapTol_; /**< soltuion gap tolerance */

private:

	double * obj_;    /**< original objective coefficients */
	double * lambda_; /**< lambda */

	CoinPackedMatrix * cpl_mat_; /**< coupling constraint matrix */
	int *    cpl_cols_;          /**< coupling columns */
	double * cpl_rhs_;           /**< right-hand sides of coupling rows */
	double   obj_offset_;        /**< constant offset in subproblem objective */

protected:

	const bool * parRelaxIntegrality_;
};

#endif /* SRC_SOLVER_DDSUB_H_ */
