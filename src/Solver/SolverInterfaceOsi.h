/*
 * SolverInterfaceOsi.h
 *
 *  Created on: Feb 9, 2015
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_SOLVERINTERFACEOSI_H_
#define SRC_SOLVER_SOLVERINTERFACEOSI_H_

#include "Solver/SolverInterface.h"

/** Coin-OR */
#include "OsiSolverInterface.hpp"

class SolverInterfaceOsi: public SolverInterface
{
public:

	/** default constructor */
	SolverInterfaceOsi(StoParam * par);

	/** copy constructor */
	SolverInterfaceOsi(SolverInterfaceOsi * si);

	/** copy constructor */
	SolverInterfaceOsi(StoParam * par, OsiSolverInterface * si);

	/** default destructor */
	virtual ~SolverInterfaceOsi();

	/** clone */
	virtual SolverInterface * clone();

	/** load problem */
	virtual void loadProblem(
			OsiSolverInterface * si,
			const char * probname = "null");

	/** load problem */
	virtual void loadProblem(
			CoinPackedMatrix * mat,
			const double * collb,
			const double * colub,
			const double * obj,
			const char * ctype,
			const double * rowlb,
			const double * rowub,
			const char * probname = "null");

	/** add row */
	virtual void addRow(int size, const int * indices, const double * vals, double lb, double ub);

	/** add cuts and return the number of cuts applied */
	virtual int addCuts(OsiCuts cuts, double effectiveness = 0.0);

	/** solve */
	virtual void solve();

	/** solution status */
	virtual STO_RTN_CODE getStatus();

	/** get number of rows */
	virtual int getNumRows() {return si_->getNumRows();}

	/** get number of columns */
	virtual int getNumCols() {return si_->getNumCols();}

	/** get column lower bounds */
	virtual const double * getColLower() {return si_->getColLower();}

	/** get column upper bounds */
	virtual const double * getColUpper() {return si_->getColUpper();}

	/** get row lower bounds */
	virtual const double * getRowLower() {return si_->getRowLower();}

	/** get row upper bounds */
	virtual const double * getRowUpper() {return si_->getRowUpper();}

	/** get global primal bound (upper bound in minimization) */
	virtual double getPrimalBound();

	/** get dual bound (lower bound in minimization) */
	virtual double getDualBound();

	/** get solution values */
	virtual const double * getSolution();

	/** get number of iterations */
	virtual int getIterationCount() {return si_->getIterationCount();}

	/** get number of nodes */
	virtual int getNumNodes() {return 0;}

	/** get print out level */
	virtual int getPrintLevel() {return si_->messageHandler()->logLevel();}

	/** get warm start information */
	virtual CoinWarmStart * getWarmStart() {return ws_;}

	/** set objective coefficients */
	virtual void setObjCoef(double * obj);

	/** set objective coefficients */
	virtual void setObjSense(int sense) {si_->setObjSense(sense);}

	/** set column bounds */
	virtual void setColLower(int index, double lb) {si_->setColLower(index, lb);}

	/** set column bounds */
	virtual void setColUpper(int index, double ub) {si_->setColUpper(index, ub);}

	/** set column bounds */
	virtual void setColBounds(int index, double lb, double ub) {si_->setColBounds(index, lb, ub);}

	/** set row bounds */
	virtual void setRowLower(int index, double lb) {si_->setRowLower(index, lb);}

	/** set row bounds */
	virtual void setRowUpper(int index, double ub) {si_->setRowUpper(index, ub);}

	/** set row bounds */
	virtual void setRowBounds(int index, double lb, double ub) {si_->setRowBounds(index, lb, ub);}

	/** set node limit */
	virtual void setNodeLimit(int limit) {}

	/** set iteration limit */
	virtual void setIterLimit(int limit);

	/** set wall time limit */
	virtual void setTimeLimit(double sec);

	/** set print out level */
	virtual void setPrintLevel(int level);

	/** set warm start information */
	virtual void setWarmStart(CoinWarmStart * ws);

	/** write MPS file */
	virtual void writeMps(const char * filename);

	/** retrieve solver interface */
	OsiSolverInterface * getOSI() {return si_;}

protected:

	/** initialize solver interface */
	virtual STO_RTN_CODE initialize();

	/** finalize solver interface */
	virtual STO_RTN_CODE finalize();

	OsiSolverInterface * si_; /**< OsiClp solver interface */
	CoinWarmStart * ws_;         /**< warmstart information */
};

#endif /* SRC_SOLVER_SOLVERINTERFACEOSI_H_ */
