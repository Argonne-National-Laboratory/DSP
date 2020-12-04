#ifndef SRC_SOLVER_DUALDECOMP_DDDROWORKERUB_H_
#define SRC_SOLVER_DUALDECOMP_DDDROWORKERUB_H_

#include "Solver/DualDecomp/DdWorkerUB.h"

/** A worker class for solving upper bounding subproblems 
 *  for the distributionally robust optimization. 
 */
class DdDroWorkerUB : public DdWorkerUB
{

	friend class DdMWSerial;
	friend class DdMWSync;
	friend class DdMWAsync;

public:
	/** A default constructor. */
	DdDroWorkerUB(
		DecModel *model, /**< model pointer */
		DspParams *par,	 /**< parameter pointer */
		DspMessage *message /**< message pointer */)
		: DdWorkerUB(model, par, message),
		  osi_dro_(NULL)
	{
	}

	/** A copy constructor. */
	DdDroWorkerUB(const DdDroWorkerUB &rhs);

	/** A default destructor. */
	virtual ~DdDroWorkerUB();

	/** A clone function */
	virtual DdDroWorkerUB *clone() const
	{
		return new DdDroWorkerUB(*this);
	}

	/** A virtual member for initializing solver. */
	virtual DSP_RTN_CODE init();

	/** A virtual member for solving problem. */
	virtual DSP_RTN_CODE solve();

	/** create problem */
	virtual DSP_RTN_CODE createProblem();

protected:
	DspOsi *osi_dro_; /**< solver interface for DRO upper bound */
};

#endif /* SRC_SOLVER_DUALDECOMP_DDDROWORKERUB_H_ */
