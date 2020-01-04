/*
 * BdWorkerMpi.h
 *
 *  Created on: Sep 19, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_BENDERS_BDWORKERMPI_H_
#define SRC_SOLVER_BENDERS_BDWORKERMPI_H_

#include "mpi.h"
/** Dsp */
#include "Solver/Benders/BdWorker.h"

class BdWorkerMpi: public BdWorker {
public:
	BdWorkerMpi(DecModel * model, DspParams * par, DspMessage * message, MPI_Comm comm);
	virtual ~BdWorkerMpi();

	/** initialize */
	int init();

	/** generate cuts */
	DSP_RTN_CODE generateCuts(int nx, int naux, const double* x, OsiCuts& cs);

protected:

	/** aggregate cuts, if necessary, based on the number of auxiliary variables
	 * introduced to the master.
	 *
	 * This may be derived to communicate cuts in MPI parallel implementation.
	 */
	virtual DSP_RTN_CODE collectCuts(int nx, int naux, double** cut, double* rhs, OsiCuts& cs);

private:

	/** Common member variables */
	MPI_Comm comm_;
	int comm_rank_;
	int comm_size_;

	int nx_;    /**< size of the solution sent from the master */
	int naux_;  /**< number of auxiliary variables in the solution */
	double* x_; /**< storage for the solution sent from the master */
	int* status_;     /**< solution status */
	int* cut_index_;  /**< subproblem index for which cut is generated */
	int* recvcounts_; /**< receive counts for statuses */
	int* displs_;     /**< displacements for statuses */
};

#endif /* SRC_SOLVER_BENDERS_BDWORKERMPI_H_ */
