/**
 * DspOsiClp.h
 *
 * 05/14/2020
 * Kibaek Kim
 */
 
#ifndef SRC_SOLVERINTERFACE_DSPOSICLP_H_
#define SRC_SOLVERINTERFACE_DSPOSICLP_H_

#include "SolverInterface/DspOsi.h"
#include "OsiClpSolverInterface.hpp"

class DspOsiClp : public DspOsi {
public:

    /** default constructor */
	DspOsiClp() {
		si_ = new OsiClpSolverInterface();
        clp_ = dynamic_cast<OsiClpSolverInterface*>(si_);
	}

	/** copy constructor */
	DspOsiClp(const DspOsiClp& rhs) {
        si_ = new OsiClpSolverInterface(*(rhs.clp_));
        clp_ = dynamic_cast<OsiClpSolverInterface*>(si_);
    }

	/** clone constructor */
	virtual DspOsiClp* clone() const {
        return new DspOsiClp(*this);
    }

    /** destructor */
	virtual ~DspOsiClp() {
        clp_ = NULL;
    }

	/** load quadratic objective */
	virtual void loadQuadraticObjective(const CoinPackedMatrix &mat) {
        clp_->getModelPtr()->loadQuadraticObjective(mat);
    }

	/** solve problem */
	virtual void solve() {
        // resolve() does not seem to recognize quadratic objective.
		si_->initialSolve();
        // si_->resolve();
	}

	virtual void use_simplex() {
		clp_->getModelPtr()->setSolveType(1);
	}

	virtual void use_barrier() {
		clp_->getModelPtr()->setSolveType(3);
	}

	/** solution statue */
	virtual int status() {
		int status = DSP_STAT_UNKNOWN;
        int status1 = clp_->getModelPtr()->status();

        if (status1 == -1) {
            status = DSP_STAT_ABORT;
        } else if (status1 == 0) {
            status = DSP_STAT_OPTIMAL;
        } else if (status1 == 1) {
            status = DSP_STAT_PRIM_INFEASIBLE;
        } else if (status1 == 2) {
            status = DSP_STAT_DUAL_INFEASIBLE;
        } else if (status1 == 3) {
            status = DSP_STAT_LIM_ITERorTIME;
        } else {
            status = DSP_STAT_ABORT;
        }
		return status;
	}

	/** set log level */
	virtual void setLogLevel(int level) {
        clp_->messageHandler()->setLogLevel(level);
		clp_->getModelPtr()->setLogLevel(level);
	}

	/** set number of cores */
	virtual void setNumCores(int num) {
        clp_->getModelPtr()->setNumberThreads(num);
	}

	/** set time limit */
	virtual void setTimeLimit(double time) {
        clp_->getModelPtr()->setMaximumSeconds(time);
    }

    OsiClpSolverInterface* clp_;    
};

#endif /* SRC_SOLVERINTERFACE_DSPOSICLP_H_ */