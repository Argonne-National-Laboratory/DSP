/*
 * TssDdMpi.h
 *
 *  Created on: Dec 10, 2014
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_TSSDDMPI_H_
#define SRC_SOLVER_TSSDDMPI_H_

/** MPI */
#include "mpi.h"

/** DSP */
#include "Solver/TssSolver.h"
#include "Solver/TssBdSub.h"
#include "Solver/TssDdMaster.h"
#include "Solver/TssDdSub.h"

//#define TSSDD_WRITE_FILE

class TssDdMpi : public TssSolver
{
	typedef vector<CoinPackedVector*> Solutions;

public:

	/** default constructor */
	TssDdMpi(
			MPI_Comm comm,
			string logfile_prefix);

	/** default destructor */
	virtual ~TssDdMpi();

	/** solve */
	virtual STO_RTN_CODE solve();

private:

	/** initialize solver */
	STO_RTN_CODE initialize();

	/** create subproblem */
	virtual STO_RTN_CODE createSubproblem();

	/** solve subproblem */
	virtual STO_RTN_CODE solveSubproblem();

	/** synchronize subproblem: send Lagrangian multipliers to subproblems */
	virtual STO_RTN_CODE syncSubproblem();

	/** create master problem */
	virtual STO_RTN_CODE createMaster();

	/** solve master problem */
	virtual STO_RTN_CODE solveMaster();

	/** synchronize master problem: send subproblem solutions to master and may update dual bound */
	virtual STO_RTN_CODE syncMaster();

	/** check whether to terminate or not */
	virtual bool doTerminate();

	/** collect results */
	virtual STO_RTN_CODE collectResults();

	/** get upper bound */
	double getUpperBound(
			TssDdSub * ddsub, /**< subproblem to evaluate */
			bool & feasible   /**< indicating feasibility */);

	/** get upper bound */
	double getUpperBound(
			const double * solution, /**< solution to evaluate */
			bool & feasible          /**< indicating feasibility */);

	/** check solution status and determine whether to continue or stop. */
	bool checkStatus(TssDdSub * ddsub);

	/** check whether solution is duplicate or not; return NULL if duplicate */
	static CoinPackedVector * duplicateSolution(
			int size,                           /**< size of array */
			const double * x,                   /**< current solution */
			vector<CoinPackedVector*> solutions /**< solution pool to check duplication */);

	/** synchronize solution pool */
	void syncSolutions(
			vector<CoinPackedVector*> & solutions,
			int & numSyncedSolutions);

	/** obtain upper bounds for solution pools */
	void obtainUpperBounds(
			int num /**< number of solutions to evaluate */);

	/** collect objective values;
	 * This takes the summation of all vals from different processors. */
	void collectObjValues(
			int num,      /**< size of vals */
			double * vals /**< objective values */);

	/** add feasibility cuts */
	int addFeasCuts(
			Solutions solutions /**< solutions to check */);

	/** add optimality cuts */
	int addOptCuts(
			int num /**< number of solutions to check */);

	/** synchronize cut pool */
	void syncCuts(
			int start,     /**< start index to synchronize */
			OsiCuts * cuts /**< cuts to synchronize */);

	/** update wall clock time remains */
	virtual void ticToc();

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

	/** solution pool */
	Solutions solutions_; /**< solutions collected from solver */

	/** cut pool */
	OsiCuts * cuts_;   /**< cut pool */
	int numSyncedCuts_; /**< number of cuts that were already synced */

	/** Upper bounding procedure */
	TssBdSub * bdsub_;
	vector<TssDdSub*> subprobs_; /**< set of subproblems (no element in root) */
	Solutions ubSolutions_;      /**< saved solutions that were evaluated for upper bounds */
	int numSyncedUbSolutions_;   /**< number of solutions that were already synced */

	/** Member variables only used in root */
	TssDdMaster * master_; /**< master problem */
	int * nsubprobs_;      /**< number of subproblems taken by each rank */
	int * scenarioSpecs_;  /**< scenario indices of which each process take care */
	double * multipliers_;  /**< Lagrangian multipliers */

	double tic_; /**< keep wall clock time for time manage */

public:

	int iterCnt_; /**< iteration counter */

	/** statistics */
	int nInfeasible_; /**< number of infeasible solutions evaluated for upper bound */
	double ctime_start_; /**< start time (cpu) */
	double wtime_start_; /**< start time (wall) */
	vector<double> wtime_elapsed_; /**< elapsed time per iteration */
	vector<double> wtime_master_;   /**< master solution time */
	vector<double> wtime_subprob_; /**< subproblem solution time */
	vector<double> objval_master_;  /**< master objective values per iteration */
	vector<double> objval_subprob_; /**< subproblem objective values per iteration */
	vector<double> primalBounds_;   /**< primal bounds per iteration */
	vector<double> dualBounds_;     /**< dual bounds per iteration */
	vector<double> changesOfMultiplier_; /**< changes of multipliers */

	string logfile_prefix_; /**< prefix of logfile name */
};

#endif /* SRC_SOLVER_TSSDDMPI_H_ */
