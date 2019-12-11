/*
 * DdWorkerLB.h
 *
 *  Created on: Apr 4, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DUALDECOMP_DDWORKERLB_H_
#define SRC_SOLVER_DUALDECOMP_DDWORKERLB_H_

#include "Solver/DualDecomp/DdWorker.h"
#include "Solver/DualDecomp/DdSub.h"

/** A worker class for solving lower bounding subproblems. */
class DdWorkerLB: public DdWorker {

	friend class DdMWSerial;
	friend class DdMWSync;
	friend class DdMWAsync;
	friend class DdMWAsyncDyn;

public:

	/** A default constructor. */
	DdWorkerLB(
			DecModel *   model,  /**< model pointer */
			DspParams *  par,    /**< parameter pointer */
			DspMessage * message /**< message pointer */);

	/** A copy constructor. */
	DdWorkerLB(const DdWorkerLB& rhs);

	/** A default destructor. */
	virtual ~DdWorkerLB();

	/** A clone function */
	virtual DdWorkerLB* clone() const {
		return new DdWorkerLB(*this);
	}

	/** A virtual member for initializing worker */
	virtual DSP_RTN_CODE init();

	/** A virtual member for solving LB problem */
	virtual DSP_RTN_CODE solve();

	/**@name Get functions */
	//@{

	/** get worker type */
	virtual int getType() {return LB;}

	/** get number of subproblems */
	virtual int getNumSubprobs() {return subprobs_.size();}

	//@}

private:

	/** A virtual member for creating LB problem */
	virtual DSP_RTN_CODE createProblem(int nsubprobs, int* subindex);

protected:

	int solution_key_; /**< solution ID to be evaluated */
	vector<DdSub*> subprobs_; /**< set of subproblems */
	bool isInit_; /**< indicate if this is the initial iteration */
};

#endif /* SRC_SOLVER_DUALDECOMP_DDWORKERLB_H_ */
