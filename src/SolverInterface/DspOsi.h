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
#include "OsiCuts.hpp"

struct callback_usr_data{
	int (*functionptr)(void *, int);
	void *cbdata;
	int *where;
};

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
	virtual void loadQuadraticObjective(const CoinPackedMatrix &mat) {
		throw CoinError("Quadratic objective is not supported.", "loadQuadraticObjective", "DspOsi");
	}

	virtual void writeMps(const char *filename){
		si_->writeMps(filename);
	}

	/** solve problem */
	virtual void solve() = 0;

	virtual void use_simplex() {
		throw CoinError("Simplex is not supported.", "use_simplex", "DspOsi");
	}

	virtual void use_barrier() {
		throw CoinError("Barrier is not supported.", "use_barrier", "DspOsi");
	}

	static int dsp_status(const OsiSolverInterface *si)
	{
		if (si->isProvenOptimal())
			return DSP_STAT_OPTIMAL;
		else if (si->isProvenPrimalInfeasible())
			return DSP_STAT_PRIM_INFEASIBLE;
		else if (si->isProvenDualInfeasible())
			return DSP_STAT_DUAL_INFEASIBLE;
		else if (si->isPrimalObjectiveLimitReached())
			return DSP_STAT_LIM_PRIM_OBJ;
		else if (si->isDualObjectiveLimitReached())
			return DSP_STAT_LIM_DUAL_OBJ;
		else if (si->isIterationLimitReached())
			return DSP_STAT_LIM_ITERorTIME;
		else if (si->isAbandoned())
			return DSP_STAT_ABORT;
		else
			return DSP_STAT_UNKNOWN;
	}

	/** solution statue */
	virtual int status() { return dsp_status(si_); }

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
	
	/** set node information display frequency */
    virtual void setNodeInfoFreq(int level) {}

	/** set number of cores */
	virtual void setNumCores(int num) {}

	/** set time limit */
	virtual void setTimeLimit(double time) {}

	/** set node limit */
	virtual void setNodeLimit(double num) {}

	/** set relative MIP gap */
	virtual void setRelMipGap(double tol) {}

	/** ==============================================================================
	  *    Generic Callbacks
	  * ============================================================================== */

	/** wrapper for cblazy */
	virtual void CallbackLazyCut(void *cbdata, OsiRowCut *lazyCut, int wherefrom=0, int purgeable=0){
		throw CoinError("Lazy callback constraints is not support", "CallbackLazyCut", "DspOsi");
	}
	virtual void CallbackLazyCut(void *cbdata, int lazylen, const int *lazyind, const double *lazyval, char lazysense, double lazyrhs){
		throw CoinError("Lazy callback constraints is not support", "CallbackLazyCut", "DspOsi");
	}
	virtual void setCallbackFunc(void (*my_callback_func)(void*)){
		throw CoinError("Lazy callback constraints is not support", "setCallbackFunc", "DspOsi");

	}
	//virtual void setCallbackFunc()
	OsiSolverInterface *si_;

	/** To store cuts added */
	OsiCuts *cutsToAdd_;
};

#endif /* SRC_SOLVERINTERFACE_DSPOSI_H_ */