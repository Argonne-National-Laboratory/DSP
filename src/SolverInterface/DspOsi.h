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

class DspOsi {
public:

	/** default constructor */
	DspOsi() :
	ismip_(false),
	isqp_(false),
	isqcp_(false) {}

	/** copy constructor */
	DspOsi(const DspOsi& rhs) : 
	ismip_(rhs.ismip_),
	isqp_(rhs.isqp_),
	isqcp_(rhs.isqcp_) {}

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

	/** load quadratic constrs */
	virtual void addQuadraticRows(int nqrows, int * linnzcnt, int * quadnzcnt, double * rhs, int * sense, int ** linind, double ** linval, int ** quadrow, int ** quadcol, double ** quadval){};
	
	/** throw error */
	virtual inline void checkDspOsiError(int err, std::string solverfuncname, std::string dsposimethod){};

	/** write problem file */
	virtual void writeProb(char const * filename_str, char const * filetype_str){};

	/** change problem type to MIQCP */
	virtual void switchToMIQCP(void){};

	/** change problem type to MIQP */
	virtual void switchToMIQP(void){};

	/** change problem type to QCP */
	virtual void switchToQCP(void){};

	/** change problem type to QP */
	virtual void switchToQP(void){};

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

	/** get ismip_ */
	virtual bool isMip() {return ismip_;}

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
	virtual void setNodeLimit(int num) {}

	/** set relative MIP gap */
	virtual void setRelMipGap(double tol) {}

	OsiSolverInterface *si_;
	
	/** Stores whether CPLEX' prob type is currently set to mixed-integer program */
	mutable bool ismip_;
	/** Stores whether CPLEX' prob type is currently set to (MI)QP or (MI)QCP */
	mutable bool isqp_;
	mutable bool isqcp_;
};

#endif /* SRC_SOLVERINTERFACE_DSPOSI_H_ */