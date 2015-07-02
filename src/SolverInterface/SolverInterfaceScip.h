/*
 * SolverInterfaceScip.h
 *
 *  Created on: Dec 8, 2014
 *      Author: kibaekkim
 */

#ifndef SOLVERINTERFACESCIP_H_
#define SOLVERINTERFACESCIP_H_

#include <vector>

/** DSP */
#include "SolverInterface/SolverInterface.h"

/** COIN */
#include "OsiCuts.hpp"

/** SCIP */
#include "objscip/objconshdlr.h"
#include "objscip/objbranchrule.h"

using namespace std;

class SolverInterfaceScip : public SolverInterface
{
public:

	/** default constructor */
	SolverInterfaceScip(StoParam * par);

	/** copy constructor */
	SolverInterfaceScip(SolverInterfaceScip * si);

	/** default destructor */
	virtual ~SolverInterfaceScip();

	/** clone */
	virtual SolverInterface * clone();

protected:

	/** initialize solver interface */
	virtual STO_RTN_CODE initialize();

	/** finalize solver interface */
	virtual STO_RTN_CODE finalize();

public:

	/**
	 * Model functions
	 */

	/** load problem */
	virtual void loadProblem(
			OsiSolverInterface * si,
			const char * probname = "scipprob");

	/** load problem */
	virtual void loadProblem(
			CoinPackedMatrix * mat,
			const double * collb,
			const double * colub,
			const double * obj,
			const char * ctype,
			const double * rowlb,
			const double * rowub,
			const char * probname = "scipprob");

	/** add row */
	virtual void addRow(int size, const int * indices, const double * vals, double lb, double ub);

	/** delete row */
	void delRow(int index);

	/** add cuts */
	virtual int addCuts(OsiCuts cuts, double effectivess)
	{
		printf("Warning: SolverInterfaceSCIP does not support addCuts();\n");
		return -1;
	}

	/** add constraint handler */
	virtual void addConstraintHandler(
			scip::ObjConshdlr * objconshdlr,
			bool deleteobject,
			bool isDual = false);

	/** find constriant handler */
	virtual scip::ObjConshdlr * findObjConshdlr(const char * name);

	/** add branch rule */
	virtual void addBranchrule(
			scip::ObjBranchrule * objbranchrule,
			bool deleteobject);

	/** find branch rule */
	virtual scip::ObjBranchrule * findObjBranchrule(const char * name);

	/** transform problem */
	virtual void transformProb();

	/** transform and presolve */
	virtual void presolve();

	/** free solve */
	void freeSolve(bool restart);

	/** free transformed problem */
	virtual void freeTransform();

#if 0
	/** enable upper bounding cuts */
	virtual void enableUpperBoundingCuts(bool yesNo, bool hasUb);
#endif

	/**
	 * Solution functions
	 */

	/** solve */
	virtual void solve();

	/** solution status */
	virtual STO_RTN_CODE getStatus();

	/**
	 * Get functions
	 */

	/** get number of rows */
	virtual int getNumRows() {return SCIPgetNOrigConss(scip_);}

	/** get number of columns */
	virtual int getNumCols() {return SCIPgetNOrigVars(scip_);}

	/** get column lower bounds */
	virtual const double * getColLower() {return clbd_;}

	/** get column upper bounds */
	virtual const double * getColUpper() {return cubd_;}

	/** get row lower bounds */
	virtual const double * getRowLower() {return rlbd_;}

	/** get row upper bounds */
	virtual const double * getRowUpper() {return rubd_;}

	/** get variable objective function */
	virtual double getObjCoeff(int j) {return SCIPvarGetObj(vars_[j]);}

	/** get global primal bound (upper bound in minimization) */
	virtual double getPrimalBound() {return SCIPgetPrimalbound(scip_);}

	/** get dual bound (lower bound in minimization) */
	virtual double getDualBound() {return SCIPgetDualbound(scip_);}

	/** get solution values */
	virtual const double * getSolution();

	/** get array of solutions; returns the number of solutions */
	virtual int getSolutions(double *** solutions);

	/** get number of iterations */
	virtual int getIterationCount() {return SCIPgetNLPIterations(scip_);}

	/** get number of nodes */
	virtual int getNumNodes() {return SCIPgetNTotalNodes(scip_);}

	/** get print out level */
	virtual int getPrintLevel();

	/**
	 * Set functions
	 */

	/** set objective coefficients */
	virtual void setObjCoef(double * obj);

	/** set objective coefficients */
	virtual void setObjSense(int sense);

	/** set print out level */
	virtual void setPrintLevel(int level);

	/** set column bounds */
	virtual void setColLower(int index, double lb);

	/** set column bounds */
	virtual void setColUpper(int index, double ub);

	/** set column bounds */
	virtual void setColBounds(int index, double lb, double ub);

	/** set row bounds */
	virtual void setRowLower(int index, double lb);

	/** set row bounds */
	virtual void setRowUpper(int index, double ub);

	/** set row bounds */
	virtual void setRowBounds(int index, double lb, double ub);

	/** set node limit */
	virtual void setNodeLimit(int limit);

	/** set iteration limit */
	virtual void setIterLimit(int limit) {}

	/** set time limit */
	virtual void setTimeLimit(double sec);

	/** set solution */
	virtual void setSolution(double * solution);

	/** set clock type */
	virtual void setClockType(int type);

	/** dual reduction */
	SCIP_RETCODE setDualReduction(bool yesNo);

	/**
	 * Cut functions
	 */

	/** get cuts */
	virtual const OsiCuts * getCuts();

	/** set cuts */
	virtual void setCuts(OsiCuts * cuts);

	/** clear cuts */
	virtual void clearCuts();

	/**
	 * SCIP specific functions
	 */

	/** get SCIP pointer */
	SCIP * getSCIP() {return scip_;}

	/** get SCIP variables */
	virtual SCIP_Var ** getSCIPvars() {return vars_;}

	/**
	 * Misc.
	 */

	/** write MPS file */
	virtual void writeMps(const char * filename);

private:

	SCIP * scip_;
	double * solution_;

	/** store original variables */
	int nvars_;
	SCIP_Var ** vars_;
	double * clbd_;
	double * cubd_;

	/** store original rows */
	int nconss_;
	SCIP_CONS ** conss_;
	double * rlbd_;
	double * rubd_;
};

#endif /* SOLVERINTERFACESCIP_H_ */
