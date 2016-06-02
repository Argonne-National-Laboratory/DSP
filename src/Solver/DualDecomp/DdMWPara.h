/*
 * DdMWPara.h
 *
 *  Created on: Apr 8, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DUALDECOMP_DDMWPARA_H_
#define SRC_SOLVER_DUALDECOMP_DDMWPARA_H_

#include "Utility/DspMpi.h"
#include <Solver/DualDecomp/DdMW.h>

class DdMWPara: public DdMW {
public:

	/** constructor */
	DdMWPara(
			MPI_Comm          comm,   /**< MPI communicator */
			DdMaster *        master, /**< master problem */
			vector<DdWorker*> worker  /**< worker for finding lower bounds */);

	virtual ~DdMWPara();

protected:

	/** initialize */
	virtual DSP_RTN_CODE init();

	/** finalize */
	virtual DSP_RTN_CODE finalize();

protected:

	/** Common member variables */
	MPI_Comm comm_;
	int comm_rank_;
	int comm_size_;
	bool sync_; /**< indicate whether parallelism is synchronous or asynchronous */

	enum {
		comm_color_main = 0, /**< color for a master and workers that solve the Lagrangian subproblems */
		comm_color_cg,
		comm_color_ub
	};

	MPI_Comm subcomm_;     /**< split communicator */
	int subcomm_size_;     /**< number of processors in each subcommunicator */
	int num_comm_colors_;  /**< number of rank colors */
	int comm_color_;       /**< rank color for subcommunicator */
	int comm_key_;         /**< rank in subcommunicator */
	int * comm_root_keys_; /**< array of ranks for each root of subcomm_ */

	MPI_Comm  lb_comm_;  /**< communicator for lower bounds */
	MPI_Group lb_group_; /**< group for lower bounds */
	MPI_Comm  cg_comm_;  /**< communicator for cut generators */
	MPI_Group cg_group_; /**< group for cut generators */
	int cg_comm_rank_;   /**< rank for cut generators */
	MPI_Comm  ub_comm_;  /**< communicator for upper bounds */
	MPI_Group ub_group_; /**< group for upper bounds */
	int ub_comm_rank_;   /**< rank for upper bounds */

	int * nsubprobs_;       /**< number of subproblems for each process */
	int * subprob_indices_; /**< subproblem indices for each process */
	int * subprob_displs_;  /**< displacement of subproblem indices */

//	bool splitWorkers_; /**< split workers for lower-bounding and upper-bounding */
};

#endif /* SRC_SOLVER_DUALDECOMP_DDMWPARA_H_ */
