/*
 * DdWorkerUB.h
 *
 *  Created on: Mar 28, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DUALDECOMP_DDWORKERUB_H_
#define SRC_SOLVER_DUALDECOMP_DDWORKERUB_H_

#include "Solver/DualDecomp/DdWorker.h"
#include "SolverInterface/DspOsi.h"

/** A worker class for solving upper bounding subproblems. */
class DdWorkerUB: public DdWorker {

	friend class DdMWSerial;
	friend class DdMWSync;
	friend class DdMWAsync;

public:

	/** A default constructor. */
	DdWorkerUB(
			DecModel *   model,  /**< model pointer */
			DspParams *  par,    /**< parameter pointer */
			DspMessage * message /**< message pointer */);

	/** A copy constructor. */
	DdWorkerUB(const DdWorkerUB& rhs);

	/** A default destructor. */
	virtual ~DdWorkerUB();

	/** A clone function */
	virtual DdWorkerUB* clone() const {
		return new DdWorkerUB(*this);
	}

	/** A virtual member for initializing solver. */
	virtual DSP_RTN_CODE init();

	/** A virtual member for solving problem. */
	virtual DSP_RTN_CODE solve();

	/** A virtual memeber for finalizing solver. */
	virtual DSP_RTN_CODE finalize() {return DSP_RTN_OK;}

public:

	/** evaluate solution */
	double evaluate(CoinPackedVector * solution);
	double evaluate(int n, double * solution);

	virtual int getType() {return UB;}

protected:

	/** create DspOsi */
	virtual DspOsi * createDspOsi();

	/** create problem */
	virtual DSP_RTN_CODE createProblem();

public:
	double bestub_; /**< best upper bound */
	std::vector<std::vector<double> > primsols_; /**< primal solution for each subproblem */

private:

	CoinPackedMatrix ** mat_mp_;
	double** rlbd_org_; /**< original row lower bounds for each subproblem */
	double** rubd_org_; /**< original row upper bounds for each subproblem */

	DspOsi ** osi_;    /**< solver interface for each subproblem */
	DspOsi * osi_dro_; /**< solver interface for DRO upper bound */
	double ub_; /**< upper bound */
};

#endif /* SRC_SOLVER_DUALDECOMP_DDWORKERUB_H_ */
