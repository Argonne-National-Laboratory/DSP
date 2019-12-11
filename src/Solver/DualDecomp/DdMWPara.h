/*
 * DdMWPara.h
 *
 *  Created on: Apr 8, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DUALDECOMP_DDMWPARA_H_
#define SRC_SOLVER_DUALDECOMP_DDMWPARA_H_

#include "Utility/DspMpi.h"
#include "DdMW.h"

/** A master-worker class for parallel dual decomposition */
class DdMWPara: public DdMW {
public:

	/** A default constructor. */
	DdMWPara(
			MPI_Comm     comm,   /**< MPI communicator */
			DecModel *   model,  /**< model pointer */
			DspParams *  par,    /**< parameters */
			DspMessage * message /**< message pointer */);

	/** A copy constructor. */
	DdMWPara(const DdMWPara& rhs);

	/** A default destructor. */
	virtual ~DdMWPara();

	/** A clone function */
	virtual DdMWPara* clone() const {
		return new DdMWPara(*this);
	}

	/** A virtual member for initializing the framework. */
	virtual DSP_RTN_CODE init();

	/** A virtual memeber for finalizing the framework. */
	virtual DSP_RTN_CODE finalize();

protected:

	/** generate Benders cuts */
	DSP_RTN_CODE generateBendersCuts(
			MPI_Comm  comm,      /**< communicator */
			int       comm_rank, /**< rank in comm */
			Solutions solutions, /**< solutions at which cuts are generated */
			OsiCuts & cuts       /**< cuts generated */);

private:

	/** create communication groups */
	DSP_RTN_CODE createGroups();

	/** set subproblem indices to the processors */
	DSP_RTN_CODE setSubproblemIndices();

	/** assign root key for each communication group */
	DSP_RTN_CODE setRootKeys();

protected:

	/** Common member variables */
	MPI_Comm comm_;
	int comm_rank_;
	int comm_size_;
	bool sync_; /**< indicate whether parallelism is synchronous or asynchronous */

	/** Master-LB group */
	MPI_Comm subcomm_;
	MPI_Group subcomm_group_;
	int subcomm_size_;
	int subcomm_rank_;

	/** LB worker group */
	MPI_Comm  lb_comm_;
	MPI_Group lb_group_;
	int lb_comm_size_;
	int lb_comm_rank_;
	int lb_comm_root_;

	/** CG-UB worker group */
	bool has_cgub_comm_;
	MPI_Comm cgub_comm_;
	MPI_Group cgub_group_;
	int cgub_comm_size_;
	int cgub_comm_rank_;
	int cgub_comm_root_;

	int * nsubprobs_;       /**< number of subproblems for each process */
	int * subprob_indices_; /**< subproblem indices for each process */
	int * subprob_displs_;  /**< displacement of subproblem indices */

//	bool splitWorkers_; /**< split workers for lower-bounding and upper-bounding */
};

#endif /* SRC_SOLVER_DUALDECOMP_DDMWPARA_H_ */
