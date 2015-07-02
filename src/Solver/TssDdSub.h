/*
 * TssDdSub.h
 *
 *  Created on: Dec 10, 2014
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_TSSDDSUB_H_
#define SRC_SOLVER_TSSDDSUB_H_

/** DSP */
#include "SolverInterface/SolverInterface.h"
#include "Solver/TssBdSub.h"
#include "Model/TssModel.h"

class TssDdSub
{
public:

	/** default constructor */
	TssDdSub(int s, StoParam * par) :
		par_(par),
		si_(NULL),
		sind_(s),
		ncols_first_(0),
		obj_(NULL),
		lambda_(NULL),
		nsols_(0),
		solutions_(NULL)
	{
		/** nothing to do */
	}

	/** default destructor */
	virtual ~TssDdSub();

	/** create problem */
	STO_RTN_CODE createProblem(
			double     probability, /**< probability */
			TssModel * model);

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

	int sind_;        /**< scenario index */
	int ncols_first_; /**< number of first-stage columns */

private:

	double * obj_;    /**< original objective coefficients */
	double * lambda_; /**< lambda */

	/** solution pool */
	int nsols_; /**< number of solutions */
	double ** solutions_; /**< solutions */

public:

	vector<double> wtime_solution_; /**< solution wall time */
};

#endif /* SRC_SOLVER_TSSDDSUB_H_ */
