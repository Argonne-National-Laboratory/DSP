/*
 * DecSolver.h
 *
 *  Created on: Sep 24, 2014
 *      Author: kibaekkim, ctjandra
 */

#ifndef DECSOLVER_H_
#define DECSOLVER_H_

#include <iostream>
#include <fstream>
#include <iomanip>
#include <memory>
#include "OsiSolverInterface.hpp"
#include "Utility/DspMessage.h"
#include "Utility/DspRtnCodes.h"
#include "Utility/StoConfig.h"
#include "Utility/DspParams.h"
#include "Model/DecModel.h"
#include "TreeSearch/DspBranchObj.h"
#include "Utility/DspMpi.h"
#include "Utility/DspUtility.h"
#include "SolverInterface/DspOsi.h"

/**
 * Abstract class for a decomposition solver implementation.
 */
class DecSolver {
public:

	/** A default constructor. */
	DecSolver(
			DecModel *   model,   /**< model pointer */
			DspParams *  par,     /**< parameter pointer */
			DspMessage * message /**< message pointer */);

	/** A copy constructor. */
	DecSolver(const DecSolver& rhs);

	/** A default destructor */
	virtual ~DecSolver();

	/** A clone function */
	virtual DecSolver* clone() const = 0;

	/** A pure virtual member for initializing solver. */
	virtual DSP_RTN_CODE init() = 0;

	/** A pure virtual member for solving problem. */
	virtual DSP_RTN_CODE solve() = 0;

	/** A pure virtual memeber for finalizing solver. */
	virtual DSP_RTN_CODE finalize() = 0;

	/**@name Get functions */
	//@{

	/** A virtual member to return the number of columns. */
	virtual int getNumCols() {return getSiPtr()->getNumCols();}

	/** A virtual member fo return the number of rows. */
	virtual int getNumRows() {return getSiPtr()->getNumRows();}

	/** A virtual member to return column lower bound. */
	virtual const double* getColLower() {return getSiPtr()->getColLower();}

	/** A virtual member to return column lower upper. */
	virtual const double* getColUpper() {return getSiPtr()->getColUpper();}

	/** A virtual member to get best primal solution */
	virtual const double * getBestPrimalSolution() {return &bestprimsol_[0];}

	/** A virtual member to get primal solution */
	virtual const double * getPrimalSolution() {return &primsol_[0];}

	/** A virtual member to get dual solution */
	virtual const double * getBestDualSolution() {return &bestdualsol_[0];}

	/** A virtual member to get dual solution */
	virtual const double * getDualSolution() {return &dualsol_[0];}

	/** A virtual member to get best primal objective */
	virtual double getBestPrimalObjective() {return bestprimobj_;}

	/** A virtual member to get primal objective */
	virtual double getPrimalObjective() {return primobj_;}

	/** A virtual member to get best dual objective */
	virtual double getBestDualObjective() {return bestdualobj_;}

	/** A virtual member to get dual objective */
	virtual double getDualObjective() {return dualobj_;}

	/** A virtual member to return absolute duality gap */
	virtual double getAbsDualityGap() {return fabs(bestprimobj_-bestdualobj_);}

	/** A virtual member to return relative duality gap */
	virtual double getRelDualityGap() {return fabs(bestprimobj_-bestdualobj_) / (1.0e-10 + fabs(bestprimobj_));}

	/** A virtual member to return absolute approximate gap */
	virtual double getAbsApproxGap() {return fabs(primobj_-bestdualobj_);}

	/** A virtual member to return relative approximate gap */
	virtual double getRelApproxGap() {return fabs(primobj_-bestdualobj_) / (1.0e-10 + fabs(primobj_));}

	/** A virtual member to get log level */
	virtual int getLogLevel() {return message_->logLevel_;}

	/** A virtual member to get model pointer */
	virtual DecModel * getModelPtr() {return model_;}

	/** A virtual member to get parameter pointer */
	virtual DspParams * getParPtr() {return par_;}

	/** A virtual member to get message pointer */
	virtual DspMessage * getMessagePtr() {return message_;}

	/** A virtual member to get solver interface */
	virtual DspOsi* getDspOsiPtr() {return osi_;}

	/** A virtual member to get solver interface */
	virtual OsiSolverInterface* getSiPtr() {return osi_->si_;}

	/** A virtual member to get solution time */
	virtual double getCpuTime() {return cputime_;}
	virtual double getWallTime() {return walltime_;}

	/** A virtual member to get number of iterations */
	virtual int getNumIterations() {return numIterations_;}

	/** A virtual member to get number of branch-and-bound nodes */
	virtual int getNumNodes() {return numNodes_;}

	/** get status */
	virtual DSP_RTN_CODE getStatus() {return status_;}

	//@}

	/**@name Set functions */
	//@{

	virtual void setBranchingObjects(const DspBranchObj* branchobj) {}

	/** set column bounds */
	virtual void setColBounds(int j, double clbd, double cubd) {getSiPtr()->setColBounds(j, clbd, cubd);}
	virtual void setColBounds(const double* clbd, const double* cubd) {
		for (int j = 0; j < getSiPtr()->getNumCols(); ++j)
			getSiPtr()->setColBounds(j, clbd[j], cubd[j]);
	}

	/** set best primal objective */
	virtual void setBestPrimalObjective(double primobj) {bestprimobj_ = primobj;}

	/** set primal objective */
	virtual void setPrimalObjective(double primobj) {primobj_=primobj;}

	/** set best primal solution */
	virtual void setBestPrimalSolution(std::vector<double>& primsol) {bestprimsol_=primsol;}

	/** set dual objective */
	virtual void setBestDualObjective(double obj) {bestdualobj_=obj;}

	/** set dual objective */
	virtual void setDualObjective(double dualobj) {dualobj_=dualobj;}

	/** set best dual solution */
	virtual void setBestDualSolution(std::vector<double>& dualsol) {bestdualsol_=dualsol;}

	/** set time limit*/
	virtual void setTimeLimit(double t) {time_remains_ = t;}

	/** set iteration limit */
	virtual void setIterLimit(int n) {iterlim_ = n;}

	/** set log level */
	virtual void setLogLevel(int level) {message_->logLevel_ = level;}

	/** set status */
	void setStatus(int status) {status_ = status;}

	//@}

	/** write output to a file */
	virtual void write(const char * filename);

protected:

	/** update time stamp and time remains */
	virtual void tic() {
		tic_ = CoinGetTimeOfDay();
	}

	/** update time stamp and time remains */
	virtual void ticToc() {
		time_remains_ -= CoinGetTimeOfDay() - tic_;
		tic_ = CoinGetTimeOfDay();
	}

protected:

	DecModel * model_;     /**< DecModel object */
	DspParams * par_;      /**< parameters */
	DspMessage * message_; /**< message */

	// OsiSolverInterface* si_; /**< Coin-Osi */
	DspOsi * osi_;

	DSP_RTN_CODE status_;  /**< solution status */
	std::vector<double> bestprimsol_; /**< best primal solution */
	std::vector<double> primsol_;     /**< primal solution */
	std::vector<double> bestdualsol_; /**< bestdual solution */
	std::vector<double> dualsol_;     /**< dual solution */
	double bestprimobj_;   /**< best primal objective */
	double primobj_;       /**< primal objective */
	double bestdualobj_;   /**< best dual objective */
	double dualobj_;       /**< dual objective */

public:

	double cputime_;    /**< cpu time */
	double walltime_;   /**< wall time */
	double time_remains_;  /**< time limit */
	double tic_;           /**< time stamp */
	int numIterations_; /**< number of iterations for a given method */
	int numNodes_;      /**< number of branch-and-bound tree nodes */
	int iterlim_; /**< iteration limits */

	/** solver statistics */
	vector<DSP_RTN_CODE> s_statuses_;  /**< history of solution statuses */
	vector<double>       s_primobjs_;  /**< history of primal objective values */
	vector<double>       s_dualobjs_;  /**< history of dual objective values */
	vector<double*>      s_primsols_;  /**< history of primal solutions */
	vector<double>       s_cputimes_;  /**< history of cpu times */
	vector<double>       s_walltimes_; /**< history of wall times */
};

#endif /* DECSOLVER_H_ */
