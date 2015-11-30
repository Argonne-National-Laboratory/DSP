/*
 * TssBdMpi.h
 *
 *  Created on: Nov 19, 2015
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_TSSBDMPI_H_
#define SRC_SOLVER_TSSBDMPI_H_

/** MPI */
#include "mpi.h"

/** DSP */
#include "Solver/TssBd.h"

class TssBdMpi: public TssBd {
public:

	enum
	{
		MASTER_NEEDS_CUTS = 0,
		MASTER_STOPPED
	};

	/** default constructor */
	TssBdMpi(MPI_Comm comm);

	/** default destructor */
	virtual ~TssBdMpi();

	/** solve */
	virtual STO_RTN_CODE solve();

	/** set auxiliary variable data */
	void setAuxColData(int size, double * obj, double * clbd, double * cubd);

	/** set auxiliary variable data */
	void setAugScenarios(int size, int * indices);

private:

	/** initialize solver */
	STO_RTN_CODE initialize();

	/** construct master problem */
	STO_RTN_CODE constructMasterProblem(TssBdSub * tssbdsub, double lowerbound);

	/** to find a lower bound by solving a set of group subproblems */
	STO_RTN_CODE findLowerBound(
			const double * probability,
			double & lowerbound);

	/** configure master */
	STO_RTN_CODE configureMaster(TssBdSub * tssbdsub);

	/** configure workers */
	STO_RTN_CODE configureWorkers(TssBdSub * tssbdsub);

	/** run master */
	STO_RTN_CODE runMaster(TssBdSub * tssbdsub);

	/** run workers */
	STO_RTN_CODE runWorkers(TssBdSub * tssbdsub);

	/** configure Stochastic LP */
	STO_RTN_CODE configureSLP();

	/** solve Stochastic LP */
	STO_RTN_CODE solveSLP(TssBdSub * tssbdsub);

	/** configure Stochastic MILP */
	STO_RTN_CODE configureSMILP(TssBdSub * tssbdsub);

	/** solve Stochastic MILP */
	STO_RTN_CODE solveSMILP();

private:

	/** Common member variables */
	MPI_Comm comm_;
	int comm_rank_;
	int comm_size_;

	/** Number of communication groups:
	 *    If comm_size_ < (number of scenarios), then this is always one.
	 *    If comm_size_ >= (number of scenarios)
	 *       AND comm_size_ % (number of scenarios) == comm_rank_ % (number of scenarios),
	 *    then this is ceil(comm_size_ / (number of scenarios)) + 1. */
	int num_comm_groups_;

	/** Group index for cut generations and upper bounds:
	 *    The subproblem solutions are distributed to different communication groups. */
	int comm_group_;

	/** parameters */
	int parProcIdxSize_;
	int * parProcIdx_;

	double * probability_; /**< probability of scenario */
};

#endif /* SRC_SOLVER_TSSBDMPI_H_ */
