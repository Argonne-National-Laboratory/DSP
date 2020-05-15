/**
 * DspOsi.h
 *
 * 12/12/2019
 * Kibaek Kim
 */

#ifndef SRC_SOLVERINTERFACE_DSPOSI_H_
#define SRC_SOLVERINTERFACE_DSPOSI_H_

#include "Utility/DspRtnCodes.h"
#include "OsiSolverInterface.hpp"
#include "OsiClpSolverInterface.hpp"
#include "OsiCbcSolverInterface.hpp"
#include "SolverInterface/DspOsiCpx.h"
#include "SolverInterface/DspOsiScip.h"

#ifdef DSP_HAS_OOQP
#include "SolverInterface/OsiOoqpSolverInterface.hpp"
#include "SolverInterface/OoqpEps.h"
#endif

class DspOsi {
public:

	/** default constructor */
	DspOsi() {}

	/** copy constructor */
	DspOsi(const DspOsi& rhs) {}

	/** clone constructor */
	virtual DspOsi* clone() const = 0;

	/** destructor */
	virtual ~DspOsi() {
		delete si_;
		si_ = NULL;
	}

	/** load quadratic objective */
	virtual void loadQuadraticObjective(const CoinPackedMatrix &mat) {}

	/** solve problem */
	virtual void solve() = 0;

	/** solution statue */
	virtual int status() {
		if (si_->isProvenOptimal())
			return DSP_STAT_OPTIMAL;
		else if (si_->isProvenPrimalInfeasible())
			return DSP_STAT_PRIM_INFEASIBLE;
		else if (si_->isProvenDualInfeasible())
			return DSP_STAT_DUAL_INFEASIBLE;
		else if (si_->isPrimalObjectiveLimitReached())
			return DSP_STAT_LIM_PRIM_OBJ;
		else if (si_->isDualObjectiveLimitReached())
			return DSP_STAT_LIM_DUAL_OBJ;
		else if (si_->isIterationLimitReached())
			return DSP_STAT_LIM_ITERorTIME;
		else if (si_->isAbandoned())
			return DSP_STAT_ABORT;
		else
			return DSP_STAT_UNKNOWN;
	}

	/** get primal objective value */
	virtual double getPrimObjValue() {return si_->getObjValue();}

	/** get dual objective value */
	virtual double getDualObjValue() {return si_->getObjValue();}

	/** get number of branch-and-bound nodes explored */
	virtual int getNumNodes() {return 0;}

	/** set log level */
	virtual void setLogLevel(int level) {
		si_->messageHandler()->setLogLevel(level);
	}

	/** set number of cores */
	virtual void setNumCores(int num) {}

	/** set time limit */
	virtual void setTimeLimit(double time) {}

	/** set node limit */
	virtual void setNodeLimit(double num) {}

	/** set relative MIP gap */
	virtual void setRelMipGap(double tol) {}

	OsiSolverInterface *si_;
};

#endif /* SRC_SOLVERINTERFACE_DSPOSI_H_ */