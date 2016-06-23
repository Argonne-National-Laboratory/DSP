/*
 * DecDdMpi.h
 *
 *  Created on: Dec 10, 2014
 *      Author: kibaekkim, ctjandra
 */

#ifndef SRC_SOLVER_DECDDMPI_H_
#define SRC_SOLVER_DECDDMPI_H_

/** MPI */
#include "mpi.h"

/** DSP */
#include "Solver/DecSolver.h"
#include "Solver/TssBdSub.h"
#include "Solver/DecDdMaster.h"
#include "Solver/DecDdSub.h"

//#define DECDD_WRITE_FILE

class DecDdMpi : public DecSolver
{
	typedef vector<CoinPackedVector*> Solutions;

public:

	/** default constructor */
	DecDdMpi(
			MPI_Comm comm,
			string logfile_prefix);

	/** default destructor */
	virtual ~DecDdMpi();

	/** solve */
	virtual DSP_RTN_CODE solve();

protected:

	/** initialize global settings for solver */
	virtual DSP_RTN_CODE initializeGlobal();

	/** create master problem */
	virtual DSP_RTN_CODE createMaster();

	/** create subproblem */
	virtual DSP_RTN_CODE createSubproblem();

	/** initialize statistics */
	virtual DSP_RTN_CODE initializeStatistics();

	/** update wall clock time remains */
	virtual double ticToc();

private:

	/** initialize local settings for solver */
	virtual DSP_RTN_CODE initializeLocal();

	/** set subproblem solvers */
	virtual DSP_RTN_CODE configureSubSolver();

	/** solve subproblem */
	virtual DSP_RTN_CODE solveSubproblem();

	/** synchronize subproblem: send Lagrangian multipliers to subproblems */
	virtual DSP_RTN_CODE syncSubproblem();

	/** solve master problem */
	virtual DSP_RTN_CODE solveMaster();

	/** synchronize master problem: send subproblem solutions to master and may update dual bound */
	virtual DSP_RTN_CODE syncMaster();

	/** check whether to terminate or not */
	virtual bool doTerminate();

	/** collect results */
	virtual DSP_RTN_CODE collectResults();

	/** get upper bound */
	double getUpperBound(
			DecDdSub * ddsub, /**< subproblem to evaluate */
			bool & feasible   /**< indicating feasibility */);

	/** get upper bound */
	double getUpperBound(
			const double * solution, /**< solution to evaluate */
			bool & feasible          /**< indicating feasibility */);

	/** check solution status and determine whether to continue or stop. */
	bool checkStatus(DecDdSub * ddsub);

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

protected:

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

	bool use_root_process_; /**< indicate if the root process solves subproblems. */

	double * lambda_; /**< Lagrangian multipliers */

private:

	/** solution pool */
	Solutions solutions_; /**< solutions collected from solver */

	/** cut pool */
	OsiCuts * cuts_;   /**< cut pool */
	int numSyncedCuts_; /**< number of cuts that were already synced */

	/** Member variables for subproblems */
	vector<DecDdSub*> subprobs_; /**< set of subproblems (no element in root) */

	/** Upper bounding procedure */
	TssBdSub * bdsub_;
	Solutions ubSolutions_;      /**< saved solutions that were evaluated for upper bounds */
	int numSyncedUbSolutions_;   /**< number of solutions that were already synced */

	/** Member variables only used in root */
	DecDdMaster * master_;  /**< master problem */
	int * nsubprobsAtRank_; /**< number of subproblems taken by each rank */
	int * subproblemSpecs_; /**< scenario indices of which each process take care */

	double tic_; /**< keep wall clock time for time manage */

public:

	int iterCnt_; /**< iteration counter */

	/** statistics */
	int nInfeasible_; /**< number of infeasible solutions evaluated for upper bound */
	vector<double> wtime_elapsed_;  /**< elapsed time per iteration */
	vector<double> wtime_master_;   /**< master solution time */
	vector<double> wtime_subprob_;  /**< subproblem solution time */
	vector<double> objval_master_;  /**< master objective values per iteration */
	vector<double> primal_subprob_; /**< subproblem primal objective values per iteration */
	vector<double> dual_subprob_;   /**< subproblem dual objective values per iteration */
	vector<double> primalBounds_;   /**< primal bounds per iteration */
	vector<double> dualBounds_;     /**< dual bounds per iteration */
	vector<double> changesOfMultiplier_; /**< changes of multipliers */

	string logfile_prefix_; /**< prefix of logfile name */

protected:

	/** parameters */
	int parMasterAlgo_;
	int parFeasCuts_;
	int parOptCuts_;
	int parEvalUb_;
	int parProcIdxSize_;
	int * parProcIdx_;
	int parNumCutsPerIter_;
	bool parLogDualVars_;
	double parStopTol_;
};

#endif /* SRC_SOLVER_DECDDMPI_H_ */
