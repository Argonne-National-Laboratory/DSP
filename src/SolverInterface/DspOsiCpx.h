/**
 * DspOsiCpx.h
 *
 * 12/12/2019
 * Kibaek Kim
 */
 
#ifndef SRC_SOLVERINTERFACE_DSPOSICPX_H_
#define SRC_SOLVERINTERFACE_DSPOSICPX_H_

#ifdef DSP_HAS_CPX

#include "cplex.h"
#include "OsiCpxSolverInterface.hpp"
#include "SolverInterface/DspOsi.h"

class DspOsiCpx : public DspOsi {
public:

	/** default constructor */
	DspOsiCpx() {
		si_ = new OsiCpxSolverInterface();
		cpx_ = dynamic_cast<OsiCpxSolverInterface*>(si_);
	}

	/** copy constructor */
	DspOsiCpx(const DspOsiCpx& rhs) {
        si_ = new OsiCpxSolverInterface(*(rhs.cpx_));
		cpx_ = dynamic_cast<OsiCpxSolverInterface*>(si_);
	}

	/** clone constructor */
	virtual DspOsiCpx* clone() const {
		return new DspOsiCpx(*this);
	}

	/** destructor */
	virtual ~DspOsiCpx() {
		cpx_ = NULL;
	}

	/** load quadratic objective */
	virtual void loadQuadraticObjective(const CoinPackedMatrix &mat) {
		if (mat.isColOrdered()) {
			for (int j = 0; j < mat.getMajorDim(); ++j) {
				for (int k = 0; k < mat.getVectorSize(j); ++k) {
					int i = mat.getIndices()[mat.getVectorStarts()[j] + k];
					double v = mat.getElements()[mat.getVectorStarts()[j] + k];
					CPXchgqpcoef(cpx_->getEnvironmentPtr(), cpx_->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL),
						i, j, v);
				}
			}
		} else {
			for (int i = 0; i < mat.getMajorDim(); ++i) {
				for (int k = 0; k < mat.getVectorSize(i); ++k) {
					int j = mat.getIndices()[mat.getVectorStarts()[i] + k];
					double v = mat.getElements()[mat.getVectorStarts()[i] + k];
					CPXchgqpcoef(cpx_->getEnvironmentPtr(), cpx_->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL),
						i, j, v);
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

	virtual void use_simplex() {
		CPXsetintparam(cpx_->getEnvironmentPtr(), CPX_PARAM_LPMETHOD, CPX_ALG_DUAL);
	}

	virtual void use_barrier() {
		CPXsetintparam(cpx_->getEnvironmentPtr(), CPX_PARAM_LPMETHOD, CPX_ALG_BARRIER);
		CPXsetintparam(cpx_->getEnvironmentPtr(), CPX_PARAM_BARCROSSALG, -1);
		CPXsetdblparam(cpx_->getEnvironmentPtr(), CPX_PARAM_BAREPCOMP, 1e-5);
	}

	/** solution statue */
	virtual int status() {
		int status = DSP_STAT_UNKNOWN;
		int probtype = CPXgetprobtype(cpx_->getEnvironmentPtr(), cpx_->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL));
		int stat = CPXgetstat(cpx_->getEnvironmentPtr(), cpx_->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL));

		if (probtype == CPXPROB_LP) {
			switch(stat) {
			case CPX_STAT_OPTIMAL:
			case CPX_STAT_OPTIMAL_INFEAS:
				status = DSP_STAT_OPTIMAL;
				break;
			case CPX_STAT_INFEASIBLE:
				status = DSP_STAT_PRIM_INFEASIBLE;
				break;
			case CPX_STAT_UNBOUNDED:
				status = DSP_STAT_DUAL_INFEASIBLE;
				break;
			case CPX_STAT_ABORT_OBJ_LIM:
			case CPX_STAT_ABORT_PRIM_OBJ_LIM:
				status = DSP_STAT_LIM_PRIM_OBJ;
				break;
			case CPX_STAT_ABORT_DUAL_OBJ_LIM:
				status = DSP_STAT_LIM_DUAL_OBJ;
				break;
			case CPX_STAT_ABORT_IT_LIM:
				status = DSP_STAT_STOPPED_ITER;
				break;
			case CPX_STAT_ABORT_TIME_LIM:
				status = DSP_STAT_STOPPED_TIME;
				break;
			case CPX_STAT_NUM_BEST:
			case CPX_STAT_FEASIBLE:
				status = DSP_STAT_FEASIBLE;
				break;
			case CPX_STAT_ABORT_USER:
				status = DSP_STAT_STOPPED_USER;
				break;
			default:
				status = DSP_STAT_UNKNOWN;
				break;
			}
		} else if (probtype == CPXPROB_MILP) {
			switch(stat) {
			case CPXMIP_OPTIMAL:
			case CPXMIP_OPTIMAL_TOL:
			case CPXMIP_OPTIMAL_INFEAS:
				status = DSP_STAT_OPTIMAL;
				break;
			case CPXMIP_INFEASIBLE:
			case CPXMIP_NODE_LIM_INFEAS:
			case CPXMIP_TIME_LIM_INFEAS:
				status = DSP_STAT_PRIM_INFEASIBLE;
				break;
			case CPXMIP_UNBOUNDED:
			case CPXMIP_INForUNBD:
				status = DSP_STAT_DUAL_INFEASIBLE;
				break;
			case CPXMIP_NODE_LIM_FEAS:
				status = DSP_STAT_STOPPED_NODE;
				break;
			case CPXMIP_TIME_LIM_FEAS:
				status = DSP_STAT_STOPPED_TIME;
				break;
			case CPXMIP_FAIL_FEAS:
			case CPXMIP_FAIL_INFEAS:
			case CPXMIP_MEM_LIM_FEAS:
			case CPXMIP_MEM_LIM_INFEAS:
			case CPXMIP_ABORT_FEAS:
			case CPXMIP_ABORT_INFEAS:
			case CPXMIP_FAIL_FEAS_NO_TREE:
			case CPXMIP_FAIL_INFEAS_NO_TREE:
				status = DSP_STAT_ABORT;
				break;
			default:
				status = DSP_STAT_UNKNOWN;
				break;
			}
		}
		return status;
	}

	/** get dual objective value */
	virtual double getDualObjValue() {
		double val;
		if (si_->getNumIntegers() > 0) {
			CPXgetbestobjval(cpx_->getEnvironmentPtr(), cpx_->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL), &val);
			val *= cpx_->getObjSense();
		} else {
			val = si_->getObjValue();
		}
		return val;
	}

	/** get number of branch-and-bound nodes explored */
	virtual int getNumNodes() {
		return CPXgetnodecnt(cpx_->getEnvironmentPtr(), cpx_->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL));
	}

	/** set number of cores */
	virtual void setNumCores(int num) {
		CPXsetintparam(cpx_->getEnvironmentPtr(), CPX_PARAM_THREADS, num);
	}

	/** set time limit */
	virtual void setTimeLimit(double time) {
		CPXsetdblparam(cpx_->getEnvironmentPtr(), CPX_PARAM_TILIM, time);
	}

	/** set node limit */
	virtual void setNodeLimit(double num) {
		CPXsetintparam(cpx_->getEnvironmentPtr(), CPX_PARAM_NODELIM, num);
	}

	/** set relative MIP gap */
	virtual void setRelMipGap(double tol) {
		CPXsetdblparam(cpx_->getEnvironmentPtr(), CPX_PARAM_EPGAP, tol);
	}

    OsiCpxSolverInterface* cpx_;   
};

#endif

#endif /* SRC_SOLVERINTERFACE_DSPOSICPX_H_ */