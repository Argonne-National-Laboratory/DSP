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
			bestprimsol_(NULL),
			primsol_(NULL),
			dualsol_(NULL),
			bestprimobj_(COIN_DBL_MAX),
			primobj_(COIN_DBL_MAX),
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

	/** default destructor */
	virtual ~DecSolver() {
		FREE_PTR(si_);
		FREE_ARRAY_PTR(bestprimsol_);
		FREE_ARRAY_PTR(primsol_);
		for (unsigned i = 0; i < s_primsols_.size(); ++i)
			FREE_ARRAY_PTR(s_primsols_[i]);
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

	/** The function chooses branching objects and returns the pointers. */
	virtual bool chooseBranchingObjects(
			DspBranch*& branchingUp, /**< [out] branching-up object */
			DspBranch*& branchingDn  /**< [out] branching-down object */) {
		return false;
	}

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
	virtual const double * getBestPrimalSolution() {return bestprimsol_;}

	/** get primal solution */
	virtual const double * getPrimalSolution() {return primsol_;}

	/** get dual solution */
	virtual const double * getDualSolution() {return dualsol_;}

	/** get best primal objective */
	virtual double getBestPrimalObjective() {return bestprimobj_;}

	/** get primal objective */
	virtual double getPrimalObjective() {return primobj_;}

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

	/** set dual objective */
	virtual void setDualObjective(double dualobj) {dualobj_=dualobj;}

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
	virtual void write(const char * filename) {
		ofstream myfile;
		myfile.open(filename);
		myfile << "Iter";
		myfile << ",Status";
		myfile << ",Prim";
		myfile << ",Dual";
		myfile << ",Cpu";
		myfile << ",Wall";
		myfile << "\n";
		for (unsigned i = 0; i < s_statuses_.size(); ++i)
		{
			myfile << i;
			myfile << "," << s_statuses_[i];
			myfile << "," << scientific << setprecision(5) << s_primobjs_[i];
			myfile << "," << scientific << setprecision(5) << s_dualobjs_[i];
			myfile << "," << fixed << setprecision(2) << s_cputimes_[i];
			myfile << "," << fixed << setprecision(2) << s_walltimes_[i];
			myfile << "\n";
		}
		myfile.close();
	}

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
	double * bestprimsol_; /**< best primal solution */
	double * primsol_;     /**< primal solution */
	double * dualsol_;     /**< dual solution */
	double bestprimobj_;   /**< best primal objective */
	double primobj_;       /**< primal objective */
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

public:

	/** solver statistics */
	vector<DSP_RTN_CODE> s_statuses_;  /**< history of solution statuses */
	vector<double>       s_primobjs_;  /**< history of primal objective values */
	vector<double>       s_dualobjs_;  /**< history of dual objective values */
	vector<double*>      s_primsols_;  /**< history of primal solutions */
	vector<double>       s_cputimes_;  /**< history of cpu times */
	vector<double>       s_walltimes_; /**< history of wall times */
};

#endif /* DECSOLVER_H_ */
