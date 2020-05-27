/**
 * DspOsiGrb.h
 *
 * 12/12/2019
 * Kibaek Kim
 */
 
#ifndef SRC_SOLVERINTERFACE_DSPOSIGRB_H_
#define SRC_SOLVERINTERFACE_DSPOSIGRB_H_

#ifdef DSP_HAS_GRB

#include "gurobi_c.h"
#include "OsiGrbSolverInterface.hpp"
#include "SolverInterface/DspOsi.h"

class DspOsiGrb : public DspOsi {
public:

	/** default constructor */
	DspOsiGrb() {
		si_ = new OsiGrbSolverInterface();
		grb_ = dynamic_cast<OsiGrbSolverInterface*>(si_);
	}

	/** copy constructor */
	DspOsiGrb(const DspOsiGrb& rhs) {
        si_ = new OsiGrbSolverInterface(*(rhs.grb_));
		grb_ = dynamic_cast<OsiGrbSolverInterface*>(si_);
	}

	/** clone constructor */
	virtual DspOsiGrb* clone() const {
		return new DspOsiGrb(*this);
	}

	/** destructor */
	virtual ~DspOsiGrb() {
		grb_ = NULL;
	}

	/** load quadratic objective */
	virtual void loadQuadraticObjective(const CoinPackedMatrix &mat) {
		if (mat.isColOrdered()) {
			for (int j = 0; j < mat.getMajorDim(); ++j) {
				for (int k = 0; k < mat.getVectorSize(j); ++k) {
					int i = mat.getIndices()[mat.getVectorStarts()[j] + k];
					double v = mat.getElements()[mat.getVectorStarts()[j] + k];
                    GUROBI_CALL("loadQuadraticObjective", GRBupdatemodel(grb_->getLpPtr(OsiGrbSolverInterface::KEEPCACHED_ALL)));
                    GUROBI_CALL("loadQuadraticObjective", GRBaddqpterms(grb_->getLpPtr(OsiGrbSolverInterface::KEEPCACHED_ALL),
                    1, i, j, v));
				}
			}
		} else {
			for (int i = 0; i < mat.getMajorDim(); ++i) {
				for (int k = 0; k < mat.getVectorSize(i); ++k) {
					int j = mat.getIndices()[mat.getVectorStarts()[i] + k];
					double v = mat.getElements()[mat.getVectorStarts()[i] + k];
					GUROBI_CALL("loadQuadraticObjective", GRBupdatemodel(grb_->getLpPtr(OsiGrbSolverInterface::KEEPCACHED_ALL)));
                    GUROBI_CALL("loadQuadraticObjective", GRBaddqpterms(grb_->getLpPtr(OsiGrbSolverInterface::KEEPCACHED_ALL),
                    1, i, j, v));
				}
			}
		}
    }

	/** solve problem */
	virtual void solve() {
		if (si_->getNumIntegers() > 0)
			si_->branchAndBound();
		else
			si_->initialSolve();
	}

	/** solution statue */
	virtual int status() {
		int status = DSP_STAT_UNKNOWN;
		int stat;
        GUROBI_CALL("status", GRBupdatemodel(grb_->getLpPtr()));
        GUROBI_CALL("status", GRBgetintattr(grb->getLpPtr(), "Status", &stat));

        switch(stat) {
        case GRB_OPTIMAL:
        case GRB_SUBOPTIMAL:
            status = DSP_STAT_OPTIMAL;
            break;
        case GRB_INFEASIBLE:
        case GRB_CUTOFF:
            status = DSP_STAT_PRIM_INFEASIBLE;
            break;
        case GRB_INF_OR_UNBD:
            status = DSP_STAT_DUAL_INFEASIBLE;
            break;
        case GRB_UNBOUNDED:
            status = DSP_STAT_DUAL_INFEASIBLE;
            break;
        case GRB_USER_OBJ_LIMIT:
            status = DSP_STAT_LIM_PRIM_OBJ;
            break;
        case GRB_ITERATION_LIMIT:
            status = DSP_STAT_STOPPED_ITER;
            break;
        case GRB_NODE_LIMIT:
            status = DSP_STAT_STOPPED_NODE;
            break;
        case GRB_TIME_LIMIT:
            status = DSP_STAT_STOPPED_TIME;
            break;
        case GRB_NUMERIC:
        case GRB_INTERRUPTED:
            status = DSP_STAT_STOPPED_USER;
            break;
        case GRB_LOADED:
        case GRB_INPROGRESS:
            status = DSP_STAT_ABORT;
            break;
	    default:
            status = DSP_STAT_UNKNOWN;
            break;
        return status;
	}

	/** get dual objective value */
	virtual double getDualObjValue() {
		double val;
        GUROBI_CALL("getDualObjVal", GRBupdatemodel(grb_->getLpPtr()));
		GUROBI_CALL("getDualObjVal", GRBgetdblattr(grb_->getLpPtr(), GRB_DBL_ATTR_OBJBOUND, &val));
		return val * grb_->getObjSense();
	}

	/** get number of branch-and-bound nodes explored */
	virtual int getNumNodes() {
		double node;
        GUROBI_CALL("getNumNodes", GRBupdatemodel(grb_->getLpPtr()));
        GUROBI_CALL("getNumNodes", GRBgetdblattr(grb->getLpPtr(), GRB_DBL_ATTR_NODECOUNT, &node));
        return (int)node;
	}

	/** set number of cores */
	virtual void setNumCores(int num) {
        GUROBI_CALL("setNumCores", GRBupdatemodel(grb_->getLpPtr()));
		GUROBI_CALL("setNumCores", GRBsetintparam(grb->getEnvironmentPtr(), GRB_INT_PAR_THREADS, num));
	}

	/** set time limit */
	virtual void setTimeLimit(double time) {
        GUROBI_CALL("setTimeLimit", GRBupdatemodel(grb_->getLpPtr()));
		GUROBI_CALL("setTimeLimit", GRBsetdblparam(grb_->getEnvironmentPtr(), GRB_DBL_PAR_TIMELIMIT, time));
	}

	/** set node limit */
	virtual void setNodeLimit(double num) {
        GUROBI_CALL("setNodeLimit", GRBupdatemodel(grb_->getLpPtr()));
		GUROBI_CALL("setNodeLimit", GRBsetdblparam(grb_->getEnvironmentPtr(), GRB_DBL_PAR_NODELIMIT, num));
	}

	/** set relative MIP gap */
	virtual void setRelMipGap(double tol) {
        GUROBI_CALL("setRelMipGap", GRBupdatemodel(grb_->getLpPtr()));
		GUROBI_CALL("setRelMipGap", GRBsetdblparam(grb_->getEnvironmentPtr(), GRB_DBL_PAR_MIPGAP, tol));
	}

    OsiGrbSolverInterface* grb_;   
};

#endif

#endif /* SRC_SOLVERINTERFACE_DSPOSICPX_H_ */