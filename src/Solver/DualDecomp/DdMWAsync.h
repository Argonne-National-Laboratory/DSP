/*
 * DdMWAsync.h
 *
 *  Created on: Mar 23, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DUALDECOMP_DDMWASYNC_H_
#define SRC_SOLVER_DUALDECOMP_DDMWASYNC_H_

#include <deque>
#include <Solver/DualDecomp/DdMWPara.h>

class DdMWAsync: public DdMWPara {

public:

	/** constructor */
	DdMWAsync(
			MPI_Comm     comm,   /**< MPI communicator */
			DecModel *   model,  /**< model pointer */
			DspParams *  par,    /**< parameters */
			DspMessage * message /**< message pointer */);

	/** destructor */
	virtual ~DdMWAsync();

	/** initialize */
	virtual DSP_RTN_CODE init();

	/** finalize */
	virtual DSP_RTN_CODE finalize();

protected:

	/** run master process */
	virtual DSP_RTN_CODE runMaster();

	/** run worker processes */
	virtual DSP_RTN_CODE runWorker();

private:

	/** run master initialization */
	DSP_RTN_CODE runMasterInit();

	/** run master core method */
	DSP_RTN_CODE runMasterCore();

	/** store coupling solution */
	DSP_RTN_CODE storeCouplingSolutions(Solutions & stored);

	/** receive coupling solutions */
	DSP_RTN_CODE recvCouplingSolutions(
			MPI_Comm comm, /**< communicator to broadcast solutions */
			int comm_rank, /**< processor rank of the given communicator */
			Solutions &solutions /**< received solution placeholder */);

	/** send master solution to workers */
	DSP_RTN_CODE sendMasterSolution(
			int solution_key,
			double * master_primsol,
			int worker_proc,
			int num_subprobs,
			int * subprobs,
			int * numCutsAdded);

	/** run worker initialization */
	DSP_RTN_CODE runWorkerInit();

	/** run worker core method */
	DSP_RTN_CODE runWorkerCore();

	/** run worker cut generation methods */
	DSP_RTN_CODE runWorkerCg(
			Solutions solutions /**< solutions at which cuts are generated */);

	/** run worker upper bound methods */
	DSP_RTN_CODE runWorkerUb(
			Solutions solutions /**< solutions to evaluate UB */);

	/** receive Benders cuts */
	DSP_RTN_CODE recvBendersCuts();

	/** receive upper bounds */
	DSP_RTN_CODE recvUpperBounds();

	/** push front queue */
	DSP_RTN_CODE pushFrontSolutionToQueue(
			double * solution /**< lambda to add */);

	/** push queue */
	DSP_RTN_CODE pushSolutionToQueue(
			double * solution /**< lambda to add */);

	/** pop queue */
	DSP_RTN_CODE popSolutionFromQueue();

	/** pop back queue */
	DSP_RTN_CODE popBackSolutionFromQueue();

private:

	enum {
		Q_NOT_ASSIGNED = 0,
		Q_ASSIGNED,
		Q_EVALUATED
	};

	int qid_counter_; /**< unique queue identification number */
	int max_queue_size_;
	std::deque<int>     q_id_;        /**< queue ID */
	std::deque<double*> q_solution_;    /**< lambdas in queue */
	std::deque<int*>    q_indicator_; /**< indicate if lambda is evaluated by a processor */
	std::deque<double>  q_objval_;    /**< objective value */

};

#endif /* SRC_SOLVER_DUALDECOMP_DDMWASYNC_H_ */
