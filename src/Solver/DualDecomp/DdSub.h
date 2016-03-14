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
#include "SolverInterface/SolverInterface.h"
#include "Solver/Benders/BdSub.h"
#include "Model/DecModel.h"

/**
 * This defines a class for solving a dual decomposition subproblem.
 */
class DdSub
{
public:

	/** default constructor */
	DdSub(int s, DspParams * par) :
		par_(par),
		si_(NULL),
		sind_(s),
		nrows_coupling_(0),
		ncols_coupling_(0),
		theta_(-COIN_DBL_MAX),
		gapTol_(0.0001),
		obj_(NULL),
		lambda_(NULL),
		nsols_(0),
		solutions_(NULL),
		cpl_mat_(NULL),
		cpl_cols_(NULL),
		cpl_rhs_(NULL),
		obj_offset_(0),
		parRelaxIntegrality_(NULL)
	{
		/** nothing to do */
	}

	/** default destructor */
	virtual ~DdSub();

	/** create problem */
	STO_RTN_CODE createProblem(DecModel * model);

	/** add cut generator */
	STO_RTN_CODE addCutGenerator(BdSub * bdsub);

	/** add branch rule */
	STO_RTN_CODE addBranchrule();

	/** change branch rule */
	STO_RTN_CODE chgBranchrule(double lb);

	/** update problem */
	STO_RTN_CODE updateProblem(
			double * lambda,
			double primal_bound = COIN_DBL_MAX);

	/** solve problem */
	STO_RTN_CODE solve();

	/** free solution process data */
	STO_RTN_CODE freeSolve(bool restart);

	/** free all solution process data */
	STO_RTN_CODE freeTransform();

	/** collect cuts */
	STO_RTN_CODE collectCuts(OsiCuts * cuts);

	/** push cuts */
	STO_RTN_CODE pushCuts(OsiCuts * cuts);

	/** collect solutions */
	STO_RTN_CODE collectSolutions();

	/** get number of solutions */
	int getNumSolutions() {return nsols_;}

	/** get solutions */
	double ** getSolutions() {return solutions_;}

	/** get primal bound of the objective value */
	double getPrimalBound() {return si_->getPrimalBound() + obj_offset_;}

	/** get dual bound of the objective value */
	double getDualBound() {return si_->getDualBound() + obj_offset_;}

	/** set solutions */
	void setSolutions(double * solution);

	/** set wall clock time limit */
	void setTimeLimit(double sec);

	/** set accuracy tolerance */
	void setGapTol(double tol);

	/** set print level */
	void setPrintLevel(int level);

	/** get MPI message buffer */
	double * MPImsgbuf();

	/** get MPI message buffer */
	STO_RTN_CODE MPImsgbuf(double * msgbuf);

public:

	DspParams * par_;
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

	/** solution pool */
	int nsols_; /**< number of solutions */
	double ** solutions_; /**< solutions */

	CoinPackedMatrix * cpl_mat_; /**< coupling constraint matrix */
	int * cpl_cols_;             /**< coupling columns */
	double * cpl_rhs_;           /**< right-hand sides of coupling rows */
	double obj_offset_;          /**< constant offset in subproblem objective */

public:

	/** solver statistics */
	vector<STO_RTN_CODE> s_statuses_;  /**< history of solution statuses */
	vector<double>       s_primobjs_;  /**< history of primal objective values */
	vector<double>       s_dualobjs_;  /**< history of dual objective values */
//	vector<double*>      s_primsols_;  /**< history of primal solutions */
	vector<double>       s_cputimes_;  /**< history of cpu times */
	vector<double>       s_walltimes_; /**< history of wall times */

protected:

	const bool * parRelaxIntegrality_;
};

#endif /* SRC_SOLVER_DDSUB_H_ */
