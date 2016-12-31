/*
 * DdWorkerUB.h
 *
 *  Created on: Mar 28, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DUALDECOMP_DDWORKERUB_H_
#define SRC_SOLVER_DUALDECOMP_DDWORKERUB_H_

#include <Solver/DualDecomp/DdWorker.h>

class DdWorkerUB: public DdWorker {

	friend class DdMWSerial;
	friend class DdMWSync;
	friend class DdMWAsync;

public:

	/** constructor */
	DdWorkerUB(DspParams * par, DecModel * model, DspMessage * message);

	/** destructor */
	virtual ~DdWorkerUB();

	/** initialize */
	virtual DSP_RTN_CODE init();

	/** solve */
	virtual DSP_RTN_CODE solve();

public:

	/** evaluate solution */
	double evaluate(CoinPackedVector * solution);

	virtual int getType() {return UB;}

protected:

	/** create problem */
	virtual DSP_RTN_CODE createProblem();

private:

	CoinPackedMatrix ** mat_mp_;
	double** rlbd_org_;
	double** rubd_org_;

	SolverInterface ** si_;
	double * objvals_;
	int * statuses_;
	double ub_; /**< upper bound */
};

#endif /* SRC_SOLVER_DUALDECOMP_DDWORKERUB_H_ */
