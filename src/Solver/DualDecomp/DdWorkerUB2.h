/*
 * DdWorkerUB2.h
 *
 *  Created on: Mar 18, 2022
 *      Author: geunyeongbyeon
 */

#ifndef SRC_SOLVER_DUALDECOMP_DDWORKERUB2_H_
#define SRC_SOLVER_DUALDECOMP_DDWORKERUB2_H_

// #include "Solver/DualDecomp/DdWorker.h"
#include "Solver/DualDecomp/DdWorkerUB.h"
#include "Solver/DualDecomp/DdSub.h"

/** A worker class for solving upper bounding subproblems. */
class DdWorkerUB2: public DdWorkerUB {

	friend class DdMWSerial;
	friend class DdMWSync;
	friend class DdMWAsync;

public:

	/** A default constructor. */
	DdWorkerUB2(
			DecModel *   model,  /**< model pointer */
			DspParams *  par,    /**< parameter pointer */
			DspMessage * message /**< message pointer */);

	/** A copy constructor. */
	DdWorkerUB2(const DdWorkerUB2& rhs);

	/** A default destructor. */
	virtual ~DdWorkerUB2();

	/** A clone function */
	virtual DdWorkerUB2* clone() const {
		return new DdWorkerUB2(*this);
	}

	/** A virtual member for initializing solver. */
	virtual DSP_RTN_CODE init();

	/** A virtual member for solving problem. */
	virtual DSP_RTN_CODE solve();

	/** A virtual memeber for finalizing solver. */
	// virtual DSP_RTN_CODE finalize() {return DSP_RTN_OK;}

	/** evaluate solution */
	virtual double evaluate(CoinPackedVector *solution);
	virtual double evaluate(int n, double *solution);

	virtual int getType() {return UB;}

protected:

	/** create problem */
	virtual DSP_RTN_CODE createProblem(int nsubprobs, int* subindex);

protected:
	vector<DdSub*> subprobs_; /**< set of subproblems */
};

#endif /* SRC_SOLVER_DUALDECOMP_DDWORKERUB2_H_ */
