/**
 * DspOsiCbc.h
 *
 * 05/14/2020
 * Kibaek Kim
 */
 
#ifndef SRC_SOLVERINTERFACE_DSPOSICBC_H_
#define SRC_SOLVERINTERFACE_DSPOSICBC_H_

#include "SolverInterface/DspOsi.h"
#include "OsiCbcSolverInterface.hpp"

class DspOsiCbc : public DspOsi {
public:

    /** default constructor */
	DspOsiCbc() : clbd_(NULL), cubd_(NULL) {
		si_ = new OsiCbcSolverInterface();
        cbc_ = dynamic_cast<OsiCbcSolverInterface*>(si_);
	}

	/** copy constructor */
	DspOsiCbc(const DspOsiCbc& rhs) {
        si_ = new OsiCbcSolverInterface(*(rhs.cbc_));
        cbc_ = dynamic_cast<OsiCbcSolverInterface*>(si_);
        clbd_ = new double [si_->getNumCols()];
        cubd_ = new double [si_->getNumCols()];
        CoinCopyN(rhs.clbd_, si_->getNumCols(), clbd_);
        CoinCopyN(rhs.cubd_, si_->getNumCols(), cubd_);
    }

	/** clone constructor */
	virtual DspOsiCbc* clone() const {
        return new DspOsiCbc(*this);
    }

    /** destructor */
	virtual ~DspOsiCbc() {
        if (clbd_) {
            delete [] clbd_;
            clbd_ = NULL;
        }
        if (cubd_) {
            delete [] cubd_;
            cubd_ = NULL;
        }
        cbc_ = NULL;
    }

	/** solve problem */
	virtual void solve() {
        cbc_->getModelPtr()->resetModel();

        if (!clbd_) {
            clbd_ = new double [si_->getNumCols()];
            CoinCopyN(si_->getColLower(), si_->getNumCols(), clbd_);
        }
        if (!cubd_) {
            cubd_ = new double [si_->getNumCols()];
            CoinCopyN(si_->getColUpper(), si_->getNumCols(), cubd_);
        }

		// si_->initialSolve();
		if (si_->getNumIntegers() > 0) {
            si_->initialSolve();
			si_->branchAndBound();
        } else
            si_->resolve();

        si_->setColLower(clbd_);
        si_->setColUpper(cubd_);
	}

	/** solution statue */
	virtual int status() {
		int status = DSP_STAT_UNKNOWN;

        int status1 = cbc_->getModelPtr()->status();
        int status2 = cbc_->getModelPtr()->secondaryStatus();

        if (status1 == -1) {
            status = DSP_STAT_ABORT;
        } else if (status1 == 0) {
            if (cbc_->isProvenOptimal() || status2 == 0)
                status = DSP_STAT_OPTIMAL;
            else if (cbc_->isProvenPrimalInfeasible() || status2 == 1)
                status = DSP_STAT_PRIM_INFEASIBLE;
            else if (cbc_->isProvenDualInfeasible() || status2 == 7)
                status = DSP_STAT_DUAL_INFEASIBLE;
        } else if (status1 == 1) {
            if (status2 == 2)
                status = DSP_STAT_STOPPED_GAP;
            else if (status2 == 3)
                status = DSP_STAT_STOPPED_NODE;
            else if (status2 == 4)
                status = DSP_STAT_STOPPED_TIME;
            else if (status2 == 5)
                status = DSP_STAT_STOPPED_USER;
            else if (status2 == 6)
                status = DSP_STAT_STOPPED_SOLUTION;
            else if (status2 == 8)
                status = DSP_STAT_STOPPED_ITER;
        } else if (status1 == 2) {
            status = DSP_STAT_ABORT;
        }
		return status;
	}

	/** get primal objective value */
	virtual double getPrimObjValue() {
        return cbc_->getModelPtr()->getObjValue();
    }

	/** get dual objective value */
	virtual double getDualObjValue() {
        return cbc_->getModelPtr()->getBestPossibleObjValue();
    }

	/** get number of branch-and-bound nodes explored */
	virtual int getNumNodes() {return cbc_->getNodeCount();}

	/** set log level */
	virtual void setLogLevel(int level) {
        cbc_->messageHandler()->setLogLevel(level);
		cbc_->getModelPtr()->setLogLevel(level);
	}

	/** set number of cores */
	virtual void setNumCores(int num) {
        cbc_->getModelPtr()->setNumberThreads(num);
	}

	/** set time limit */
	virtual void setTimeLimit(double time) {
        cbc_->getModelPtr()->setMaximumSeconds(time);
    }

	/** set node limit */
	virtual void setNodeLimit(double num) {
        cbc_->getModelPtr()->setMaximumNodes(num);
    }

	/** set relative MIP gap */
	virtual void setRelMipGap(double tol) {
        cbc_->getModelPtr()->setAllowableFractionGap(tol);
    }

    OsiCbcSolverInterface* cbc_;    

protected:

    double * clbd_;
    double * cubd_;
};

#endif /* SRC_SOLVERINTERFACE_DSPOSICBC_H_ */