/*
 * DecDdSub.h
 *
 *  Created on: Dec 10, 2014
 *      Author: kibaekkim, ctjandra
 */

#ifndef SRC_SOLVER_DECDDSUB_H_
#define SRC_SOLVER_DECDDSUB_H_

#include <typeinfo>

/** DSP */
#include "SolverInterface/SolverInterface.h"
#include "Solver/TssBdSub.h"
#include "Model/DecModel.h"

class DecDdSub
{
public:

	/** default constructor */
	DecDdSub(int s, StoParam * par) :
		par_(par),
		si_(NULL),
		sind_(s),
		nrows_coupling_(0),
		ncols_coupling_(0),
		obj_(NULL),
		lambda_(NULL),
		nsols_(0),
		solutions_(NULL),
		cpl_mat_(NULL),
		cpl_cols_(NULL)
	{
		/** nothing to do */
	}

	/** default destructor */
	virtual ~DecDdSub();

	/** create problem */
	STO_RTN_CODE createProblem(DecModel * model);

	/** add cut generator */
	STO_RTN_CODE addCutGenerator(TssBdSub * tss);

	/** change cut generator */
	STO_RTN_CODE chgCutGenerator(TssBdSub * tss);

	/** add branch rule */
	STO_RTN_CODE addBranchrule();

	/** change branch rule */
	STO_RTN_CODE chgBranchrule(double lb);

	/** update problem */
	STO_RTN_CODE updateProblem(
			double * lambda,
			double primal_bound);

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

	/** set solutions */
	void setSolutions(double * solution);

	/** set wall clock time limit */
	void setTimeLimit(double sec);

	/** set print level */
	void setPrintLevel(int level);

	/** get MPI message buffer */
	double * MPImsgbuf();

	/** get MPI message buffer */
	STO_RTN_CODE MPImsgbuf(double * msgbuf);

public:

	StoParam * par_;
	SolverInterface * si_; /**< solver interface */

	int sind_;           /**< scenario index */
	int nrows_coupling_; /**< number of coupling constraints for the subproblem (dimension of lambda) */
	int ncols_coupling_; /**< number of coupling variables for the subproblem
						   *  (dimension of the part of the solution returned to the master problem) */

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

	vector<double> wtime_solution_; /**< solution wall time */
};

#endif /* SRC_SOLVER_DECDDSUB_H_ */
