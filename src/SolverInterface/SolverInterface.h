/*
 * SolverInterface.h
 *
 *  Created on: Dec 9, 2014
 *      Author: kibaekkim
 */

#ifndef SOLVERINTERFACE_H_
#define SOLVERINTERFACE_H_

/** Coin-OR */
#include "CoinPackedMatrix.hpp"
#include "OsiCuts.hpp"

/** DSP */
#include "Utility/StoRtnCodes.h"
#include "Solver/StoParam.h"

class SolverInterface
{
public:

	/** default constructor */
	SolverInterface(StoParam * par) : par_(par) {}

	/** default destructor */
	virtual ~SolverInterface() {par_=NULL;}

	/** clone */
	virtual SolverInterface * clone() = 0;

protected:

	/** initialize solver interface */
	virtual STO_RTN_CODE initialize() = 0;

	/** finalize solver interface */
	virtual STO_RTN_CODE finalize() = 0;

public:

	/**
	 * Model functions
	 */

	/** load problem */
	virtual void loadProblem(
			CoinPackedMatrix * mat,
			const double * collb,
			const double * colub,
			const double * obj,
			const char * ctype,
			const double * rowlb,
			const double * rowub,
			const char * probname = "null") = 0;

	/** load problem */
	virtual void loadProblem(OsiSolverInterface * si,
			const char * probname = "null") = 0;

	/** add row */
	virtual void addRow(int size, const int * indices, const double * vals, double lb, double ub) = 0;

	/** add cuts and return the number of cuts applied */
	virtual int addCuts(OsiCuts cuts, double effectiveness = 0.0) = 0;

	/**
	 * Solution functions
	 */

	/** solve */
	virtual void solve() = 0;

	/** solution status */
	virtual STO_RTN_CODE getStatus() = 0;

	/**
	 * Get functions
	 */

	/** get number of rows */
	virtual int getNumRows() = 0;

	/** get number of columns */
	virtual int getNumCols() = 0;

	/** get column lower bounds */
	virtual const double * getColLower() = 0;

	/** get column upper bounds */
	virtual const double * getColUpper() = 0;

	/** get row lower bounds */
	virtual const double * getRowLower() = 0;

	/** get row upper bounds */
	virtual const double * getRowUpper() = 0;

	/** get global primal bound (upper bound in minimization) */
	virtual double getPrimalBound() = 0;

	/** get dual bound (lower bound in minimization) */
	virtual double getDualBound() = 0;

	/** get solution values */
	virtual const double * getSolution() = 0;

	/** get number of iterations */
	virtual int getIterationCount() = 0;

	/** get number of nodes */
	virtual int getNumNodes() = 0;

	/** get print out level */
	virtual int getPrintLevel() = 0;

	/**
	 * Set functions
	 */

	/** set objective coefficients */
	virtual void setObjCoef(double * obj) = 0;

	/** set objective coefficients */
	virtual void setObjSense(int sense) = 0;

	/** set column bounds */
	virtual void setColLower(int index, double lb) = 0;

	/** set column bounds */
	virtual void setColUpper(int index, double ub) = 0;

	/** set column bounds */
	virtual void setColBounds(int index, double lb, double ub) = 0;

	/** set row bounds */
	virtual void setRowLower(int index, double lb) = 0;

	/** set row bounds */
	virtual void setRowUpper(int index, double ub) = 0;

	/** set row bounds */
	virtual void setRowBounds(int index, double lb, double ub) = 0;

	/** set node limit */
	virtual void setNodeLimit(int limit) = 0;

	/** set iteration limit */
	virtual void setIterLimit(int limit) = 0;

	/** set wall time limit */
	virtual void setTimeLimit(double sec) = 0;

	/** set print out level */
	virtual void setPrintLevel(int level) = 0;

	/**
	 * Misc.
	 */

	/** write MPS file */
	virtual void writeMps(const char * filename) = 0;

public:

	StoParam * par_;
};

#endif /* SOLVERINTERFACE_H_ */
