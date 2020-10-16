/**
 * DspOsiScip.h
 *
 * 12/12/2019
 * Kibaek Kim
 */
 
#ifndef SRC_SOLVERINTERFACE_DSPOSISCIP_H_
#define SRC_SOLVERINTERFACE_DSPOSISCIP_H_

#ifdef DSP_HAS_SCIP

#include "SolverInterface/DspOsi.h"
#include "SolverInterface/OsiScipSolverInterface.hpp"
#include "SolverInterface/SCIPbranchruleLB.h"
#include "CoinHelperFunctions.hpp"

class DspOsiScip : public DspOsi {
public:

    /** default constructor */
    DspOsiScip() {
        si_ = new OsiScipSolverInterface();
        scip_ = dynamic_cast<OsiScipSolverInterface*>(si_);
    }

    /** copy constructor */
    DspOsiScip(const DspOsiScip& rhs) {
        si_ = new OsiScipSolverInterface(*(rhs.scip_));
        scip_ = dynamic_cast<OsiScipSolverInterface*>(si_);
    }

    /** clone constructor */
    virtual DspOsiScip* clone() const {
        return new DspOsiScip(*this);
    }

    /** destructor */
    virtual ~DspOsiScip() {
        scip_ = NULL;
    }

    /** solve problem */
    virtual void solve() {
        si_->branchAndBound();
    }

	virtual void use_simplex() {}

    /** solution statue */
    virtual int status() {
        int status = DSP_STAT_UNKNOWN;    
        int scipstat = SCIPgetStatus(scip_->getScip());
        switch(scipstat) {
        case SCIP_STATUS_USERINTERRUPT:
            status = DSP_STAT_STOPPED_USER;
            break;
        case SCIP_STATUS_NODELIMIT:
        case SCIP_STATUS_TOTALNODELIMIT:
            status = DSP_STAT_STOPPED_NODE;
            break;
        case SCIP_STATUS_TIMELIMIT:
            status = DSP_STAT_STOPPED_TIME;
            break;
        case SCIP_STATUS_GAPLIMIT:
            status = DSP_STAT_STOPPED_GAP;
            break;
        case SCIP_STATUS_SOLLIMIT:
        case SCIP_STATUS_BESTSOLLIMIT:
            status = DSP_STAT_STOPPED_SOLUTION;
            break;
        case SCIP_STATUS_STALLNODELIMIT:
        case SCIP_STATUS_MEMLIMIT:
        case SCIP_STATUS_RESTARTLIMIT:
            status = DSP_STAT_STOPPED_UNKNOWN;
            break;
        case SCIP_STATUS_OPTIMAL:
            status = DSP_STAT_OPTIMAL;
            break;
        case SCIP_STATUS_INFEASIBLE:
            status = DSP_STAT_PRIM_INFEASIBLE;
            break;
        case SCIP_STATUS_UNBOUNDED:
        case SCIP_STATUS_INFORUNBD:
            status = DSP_STAT_DUAL_INFEASIBLE;
            break;
        case SCIP_STATUS_UNKNOWN:
        default:
            status = DSP_STAT_UNKNOWN;
            break;
        }
        return status;
    }

    /** get dual objective value */
    virtual double getDualObjValue() {
        return SCIPgetDualbound(scip_->getScip());
    }

    /** get number of branch-and-bound nodes explored */
    virtual int getNumNodes() {return SCIPgetNNodes(scip_->getScip());}

	/** set log level */
	virtual void setLogLevel(int level) {
		SCIPsetIntParam(scip_->getScip(), "display/verblevel", CoinMax(0,CoinMin(5,level)));
	}

    /** set node information display frequency */
    virtual void setNodeInfoFreq(int level) {
        SCIPsetIntParam(scip_->getScip(), "display/freq", CoinMax(-1,level));
    }

    /** set number of cores */
    virtual void setNumCores(int num)
    {
        SCIPsetIntParam(scip_->getScip(), "parallel/maxnthreads", num);
    }

    /** set time limit */
    virtual void setTimeLimit(double time) {
        SCIPsetRealParam(scip_->getScip(), "limits/time", CoinMin(time,1.0e+20));
    }

    /** set node limit */
    virtual void setNodeLimit(double num) {
        SCIPsetLongintParam(scip_->getScip(), "limits/nodes", num);
    }

    /** set relative MIP gap */
    virtual void setRelMipGap(double tol) {
        SCIPsetRealParam(scip_->getScip(), "limits/gap", tol);
    }

    OsiScipSolverInterface* scip_;  
};

#endif /* DSP_HAS_SCIP */

#endif /* SRC_SOLVERINTERFACE_DSPOSISCIP_H_ */
