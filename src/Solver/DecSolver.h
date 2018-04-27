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
/** Coin */
#include "OsiSolverInterface.hpp"
/** Dsp */
#include "Utility/DspMessage.h"
#include "Utility/DspMpi.h"
#include "Utility/DspRtnCodes.h"
#include "Utility/StoConfig.h"
#include "Utility/DspParams.h"
#include "Model/DecModel.h"
#include "TreeSearch/DspBranch.h"

class DecSolver {
public:

	/** constructor */
	DecSolver(DecModel * model, DspParams * par, DspMessage * message) :
			model_(model),
			par_(par),
			message_(message),
			si_(NULL),
			status_(DSP_STAT_UNKNOWN),
			bestprimobj_(COIN_DBL_MAX),
			primobj_(COIN_DBL_MAX),
			bestdualobj_(-COIN_DBL_MAX),
			dualobj_(-COIN_DBL_MAX),
			absgap_(COIN_DBL_MAX),
			relgap_(COIN_DBL_MAX),
			cputime_(0.0),
			walltime_(0.0),
			time_remains_(COIN_DBL_MAX),
			tic_(0.0),
			numIterations_(0),
			numNodes_(0),
			iterlim_(COIN_INT_MAX) {
		/** nothing to do */
	}

	/** copy constructor */
	DecSolver(const DecSolver& rhs):
		model_(rhs.model_),
		par_(rhs.par_),
		message_(rhs.message_),
		si_(rhs.si_->clone()),
		status_(rhs.status_),
		bestprimobj_(rhs.bestprimobj_),
		primobj_(rhs.primobj_),
		bestdualobj_(rhs.bestdualobj_),
		dualobj_(rhs.dualobj_),
		bestprimsol_(rhs.bestprimsol_),
		primsol_(rhs.primsol_),
		bestdualsol_(rhs.bestdualsol_),
		dualsol_(rhs.dualsol_),
		absgap_(rhs.absgap_),
		relgap_(rhs.relgap_),
		cputime_(rhs.cputime_),
		walltime_(rhs.walltime_),
		time_remains_(rhs.time_remains_),
		tic_(rhs.tic_),
		numIterations_(rhs.numIterations_),
		numNodes_(rhs.numNodes_),
		iterlim_(rhs.iterlim_) {
		/** nothing to do */
	}

	/** copy operator */
	DecSolver& operator=(const DecSolver& rhs) {
		model_ = rhs.model_;
		par_ = rhs.par_;
		message_ = rhs.message_;
		si_ = rhs.si_->clone();
		status_ = rhs.status_;
		bestprimobj_ = rhs.bestprimobj_;
		primobj_ = rhs.primobj_;
		bestdualobj_ = rhs.bestdualobj_;
		dualobj_ = rhs.dualobj_;
		bestprimsol_ = rhs.bestprimsol_;
		primsol_ = rhs.primsol_;
		bestdualsol_ = rhs.bestdualsol_;
		dualsol_ = rhs.dualsol_;
		absgap_ = rhs.absgap_;
		relgap_ = rhs.relgap_;
		cputime_ = rhs.cputime_;
		walltime_ = rhs.walltime_;
		time_remains_ = rhs.time_remains_;
		tic_ = rhs.tic_;
		numIterations_ = rhs.numIterations_;
		numNodes_ = rhs.numNodes_;
		iterlim_ = rhs.iterlim_;
		return *this;
	}

	virtual DecSolver* clone() const = 0;

	/** default destructor */
	virtual ~DecSolver() {
		FREE_PTR(si_);
		message_ = NULL;
		par_ = NULL;
		model_ = NULL;
	}

	/** initialize */
	virtual DSP_RTN_CODE init() = 0;

	/** solve */
	virtual DSP_RTN_CODE solve() = 0;

	/** finalize */
	virtual DSP_RTN_CODE finalize() = 0;

public:

	/**@name Get functions */
	//@{

	/** number of columns */
	virtual int getNumCols() {return si_->getNumCols();}

	/** number of rows */
	virtual int getNumRows() {return si_->getNumRows();}

	/** column lower bound */
	virtual const double* getColLower() {return si_->getColLower();}

	/** column lower upper */
	virtual const double* getColUpper() {return si_->getColUpper();}

	/** is column j integer? */
	virtual bool isInteger(int j) {return si_->isInteger(j);}

	/** solver status */
	virtual DSP_RTN_CODE getStatus() {return status_;}

	/** get best primal solution */
	virtual const double * getBestPrimalSolution() {return &bestprimsol_[0];}

	/** get primal solution */
	virtual const double * getPrimalSolution() {return &primsol_[0];}

	/** get dual solution */
	virtual const double * getBestDualSolution() {return &bestdualsol_[0];}

	/** get dual solution */
	virtual const double * getDualSolution() {return &dualsol_[0];}

	/** get best primal objective */
	virtual double getBestPrimalObjective() {return bestprimobj_;}

	/** get primal objective */
	virtual double getPrimalObjective() {return primobj_;}

	/** get best dual objective */
	virtual double getBestDualObjective() {return bestdualobj_;}

	/** get dual objective */
	virtual double getDualObjective() {return dualobj_;}

	/** get log level */
	virtual int getLogLevel() {return message_->logLevel_;}

	/** get model pointer */
	virtual DecModel * getModelPtr() {return model_;}

	/** get parameter pointer */
	virtual DspParams * getParPtr() {return par_;}

	/** get message pointer */
	virtual DspMessage * getMessagePtr() {return message_;}

	/** get solver interface */
	virtual OsiSolverInterface* getSiPtr() {return si_;}

	/** get solution time */
	virtual double getCpuTime() {return cputime_;}
	virtual double getWallTime() {return walltime_;}

	/** get number of iterations */
	virtual int getNumIterations() {return numIterations_;}

	/** get number of branch-and-bound nodes */
	virtual int getNumNodes() {return numNodes_;}

	//@}

	/**@name Set functions */
	//@{

	virtual void setBranchingObjects(const DspBranch* branchobj) {}

	/** set column bounds */
	virtual void setColBounds(int j, double clbd, double cubd) {si_->setColBounds(j, clbd, cubd);}
	virtual void setColBounds(const double* clbd, const double* cubd) {
		for (int j = 0; j < si_->getNumCols(); ++j)
			si_->setColBounds(j, clbd[j], cubd[j]);
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

	OsiSolverInterface* si_; /**< Coin-Osi */

	DSP_RTN_CODE status_;  /**< solution status */
	std::vector<double> bestprimsol_; /**< best primal solution */
	std::vector<double> primsol_;     /**< primal solution */
	std::vector<double> bestdualsol_; /**< bestdual solution */
	std::vector<double> dualsol_;     /**< dual solution */
	double bestprimobj_;   /**< best primal objective */
	double primobj_;       /**< primal objective */
	double bestdualobj_;   /**< best dual objective */
	double dualobj_;       /**< dual objective */
	double absgap_;        /**< absolute primal-dual gap */
	double relgap_;        /**< relative primal-dual gap */

	double cputime_;    /**< cpu time */
	double walltime_;   /**< wall time */
	double time_remains_;  /**< time limit */
	double tic_;           /**< time stamp */
	int numIterations_; /**< number of iterations for a given method */
	int numNodes_;      /**< number of branch-and-bound tree nodes */

	int iterlim_; /**< iteration limits */
};

#endif /* DECSOLVER_H_ */
