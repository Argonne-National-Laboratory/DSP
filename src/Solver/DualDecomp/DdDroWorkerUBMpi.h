#ifndef SRC_SOLVER_DUALDECOMP_DDDROWORKERUBMPI_H_
#define SRC_SOLVER_DUALDECOMP_DDDROWORKERUBMPI_H_

#include "Utility/DspMpi.h"
#include "Solver/DualDecomp/DdDroWorkerUB.h"

/** A worker class for solving upper bounding subproblems 
 *  for the distributionally robust optimization
 *  in parallel. 
 */
class DdDroWorkerUBMpi : public DdDroWorkerUB
{

	friend class DdMWSerial;
	friend class DdMWSync;
	friend class DdMWAsync;

public:
	/** A default constructor. */
	DdDroWorkerUBMpi(
		MPI_Comm comm,	 /**< MPI communicator */
		DecModel *model, /**< model pointer */
		DspParams *par,	 /**< parameter pointer */
		DspMessage *message /**< message pointer */);

	/** A copy constructor. */
	DdDroWorkerUBMpi(const DdDroWorkerUBMpi &rhs)
		: DdDroWorkerUB(rhs)
	{
		comm_ = rhs.comm_;
		comm_size_ = rhs.comm_size_;
		comm_rank_ = rhs.comm_rank_;
	}

	/** A default destructor. */
	virtual ~DdDroWorkerUBMpi()
	{
		comm_ = MPI_COMM_NULL;
		comm_size_ = 0;
		comm_rank_ = -1;
	}

	/** A clone function */
	virtual DdDroWorkerUBMpi *clone() const
	{
		return new DdDroWorkerUBMpi(*this);
	}

	/** A virtual member for initializing solver. */
	virtual DSP_RTN_CODE init();

	/** A virtual member for solving problem. */
	virtual DSP_RTN_CODE solve();

private:
	MPI_Comm comm_;
	int comm_rank_;
	int comm_size_;
};

#endif /* SRC_SOLVER_DUALDECOMP_DDDROWORKERUBMPI_H_ */
