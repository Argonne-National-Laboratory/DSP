/*
 * DwSolverSerial.h
 *
 *  Created on: Dec 5, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DANTZIGWOLFE_DWSOLVERSERIAL_H_
#define SRC_SOLVER_DANTZIGWOLFE_DWSOLVERSERIAL_H_

#include <DecSolver.h>
#include <DantzigWolfe/DwMaster.h>

class DwSolverSerial: public DecSolver {
public:

    /** default constructor */
	DwSolverSerial(
			DecModel *   model,  /**< model pointer */
			DspParams *  par,    /**< parameters */
			DspMessage * message /**< message pointer */);

	virtual ~DwSolverSerial();

	/** initialize */
	virtual DSP_RTN_CODE init();

	/** solve */
	virtual DSP_RTN_CODE solve();

	/** finalize */
	virtual DSP_RTN_CODE finalize();

	/** The function chooses branching objects and returns the pointers. */
	virtual bool chooseBranchingObjects(
			DspBranch*& branchingUp, /**< [out] branching-up object */
			DspBranch*& branchingDn  /**< [out] branching-down object */);

public:

	/**
	 * Redefine get functions
	 */

	/** get number of original columns */
	virtual int getNumCols() {
		return master_->ncols_orig_;
	}

	/** get number of original rows */
	virtual int getNumRows() {
		return master_->nrows_orig_;
	}

	virtual DSP_RTN_CODE getStatus() {
		return master_->getStatus();
	}

	/** get best primal solution */
	virtual const double * getBestPrimalSolution() {
		return master_->getBestPrimalSolution();
	}

	/** get primal solution */
	virtual const double * getPrimalSolution() {
		return master_->getPrimalSolution();
	}

	/** get dual solution */
	virtual const double * getDualSolution() {
		return master_->getDualSolution();
	}

	/** get best primal objective */
	virtual double getBestPrimalObjective() {
		return master_->getBestPrimalObjective();
	}

	/** get primal objective */
	virtual double getPrimalObjective() {
		return master_->getPrimalObjective();
	}

	/** get dual objective */
	virtual double getDualObjective() {
		return master_->getDualObjective();
	}

	/** get heuristic runs on/off */
	virtual bool getHeuristicRuns() {
		return master_->getHeuristicRuns();
	}

	/** get log level */
	virtual int getLogLevel() {
		return master_->getLogLevel();
	}

	/**
	 * Redefine set functions
	 */

	/** set branching objects */
	virtual void setBranchingObjects(const DspBranch* branchobj);

	/** set best primal objective */
	virtual void setBestPrimalObjective(double primobj) {
		master_->setBestPrimalObjective(primobj);
	}

	/** set primal objective */
	virtual void setPrimalObjective(double primobj) {
		master_->setPrimalObjective(primobj);
	}

	/** set dual objective */
	virtual void setDualObjective(double dualobj) {
		master_->setDualObjective(dualobj);
	}

	/** set time limit */
	virtual void setTimeLimit(double t) {
		master_->setTimeLimit(t);
	}

	/** set iteration limit */
	virtual void setIterLimit(int n) {
		master_->setIterLimit(n);
	}

	/** set heuristic runs on/off */
	virtual void setHeuristicRuns(bool on) {
		master_->setHeuristicRuns(on);
	}

	/** set log level */
	virtual void setLogLevel(int level) {
		master_->setLogLevel(level);
	}

protected:

	DwMaster* master_;
	DwWorker* worker_;

};

#endif /* SRC_SOLVER_DANTZIGWOLFE_DWSOLVERSERIAL_H_ */
