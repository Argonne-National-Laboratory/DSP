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

	/** run master initialization */
	virtual DSP_RTN_CODE runMasterInit();

	/** run master core method */
	virtual DSP_RTN_CODE runMasterCore();

	/** choose queue element for evaluating dual variables */
	virtual bool chooseQueueElement(int& qid, double*& qsol, int& nsubprobs, int*& subindex);

	/** store coupling solution */
	virtual DSP_RTN_CODE storeCouplingSolutions(Solutions & stored);

	/** receive coupling solutions */
//	DSP_RTN_CODE recvCouplingSolutions(
//			MPI_Comm comm, /**< communicator to broadcast solutions */
//			int comm_rank, /**< processor rank of the given communicator */
//			Solutions &solutions /**< received solution placeholder */);

	/** send master solution to workers */
	virtual DSP_RTN_CODE sendMasterSolution(
			int solution_key,
			double * master_primsol,
			int worker_proc,
			int num_subprobs,
			int * subprobs,
			int * numCutsAdded);

	/** run worker initialization */
	virtual DSP_RTN_CODE runWorkerInit();

	/** run worker core method */
	virtual DSP_RTN_CODE runWorkerCore();

	/** set lower bounding workers */
	virtual DSP_RTN_CODE setWorkerLb(DdWorkerLB* workerlb, int nsubprobs, int* subindex, double* buf, double bestprimobj);

	/** run worker cut generation methods */
	virtual DSP_RTN_CODE runWorkerCg(
			Solutions solutions /**< solutions at which cuts are generated */);

	/** run worker upper bound methods */
	virtual DSP_RTN_CODE runWorkerUb(
			Solutions solutions /**< solutions to evaluate UB */);

	/** receive Benders cuts */
	virtual DSP_RTN_CODE recvBendersCuts();

	/** receive upper bounds */
	virtual DSP_RTN_CODE recvUpperBounds();

	/** push front queue */
	virtual DSP_RTN_CODE pushFrontSolutionToQueue(
			double * solution /**< lambda to add */);

	/** push queue */
	virtual DSP_RTN_CODE pushSolutionToQueue(
			double * solution /**< lambda to add */);

	/** pop queue */
	virtual DSP_RTN_CODE popSolutionFromQueue();

	/** pop back queue */
	virtual DSP_RTN_CODE popBackSolutionFromQueue();

protected:

	enum {
		Q_NOT_ASSIGNED = 0,
		Q_ASSIGNED,
		Q_EVALUATED
	};

	int qid_counter_;    /**< unique queue identification number */
	int max_queue_size_; /**< maximum queue size */

	std::deque<int>     q_id_;        /**< queue ID */
	std::deque<double*> q_solution_;  /**< lambdas in queue */
	std::deque<int*>    q_indicator_; /**< indicate if lambda is evaluated for each subproblem */
	std::deque<double>  q_objval_;    /**< objective value */

};

#endif /* SRC_SOLVER_DUALDECOMP_DDMWASYNC_H_ */
