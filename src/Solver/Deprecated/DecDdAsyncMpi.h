/*
 * DecDdAsyncMpi.h
 *
 *  Created on: Feb 4, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DECDDASYNCMPI_H_
#define SRC_SOLVER_DECDDASYNCMPI_H_

#include "Solver/DecDdMpi.h"

class DecDdAsyncMpi: public DecDdMpi {
public:

	/** default constructor */
	DecDdAsyncMpi(MPI_Comm comm);

	/** default destructor */
	virtual ~DecDdAsyncMpi();

	/** solve */
	virtual DSP_RTN_CODE solve();

protected:

	/** create master problem */
	virtual DSP_RTN_CODE createMaster();

	/** initialize local settings */
	virtual DSP_RTN_CODE initializeLocal();

private:

	/** run master process */
	virtual DSP_RTN_CODE runMaster();

	/** run worker processes */
	virtual DSP_RTN_CODE runWorkers();

	/** communicate with master process (significant only at worker processes) */
	virtual DSP_RTN_CODE commWithMaster(int & doContinue);

	/** update subproblems */
	virtual DSP_RTN_CODE updateSubproblems(int & doContinue);

	/** solve subproblems */
	virtual DSP_RTN_CODE solveSubproblems(int & doContinue);

private:

	int scount_; /**< sending message size */
	int rcount_; /**< receiving message size */
	double * sendbuf_; /**< sending message buffer */
	double * recvbuf_; /**< receiving message buffer */

	/** subproblem solution message (significant only for master) */
	int maxnumsubprobs_;    /**< maximum number of subproblems among workers */
	int * subindex_;        /**< array of subproblem indices */
	double * subprimobj_;   /**< subproblem primal objective values */
	double * subdualobj_;   /**< subproblem dual objective values */
	double ** subsolution_; /**< subproblem solution */
};

#endif /* SRC_SOLVER_DECDDASYNCMPI_H_ */
