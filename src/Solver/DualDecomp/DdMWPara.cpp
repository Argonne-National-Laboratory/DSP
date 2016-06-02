/*
 * DdMWPara.cpp
 *
 *  Created on: Apr 8, 2016
 *      Author: kibaekkim
 */

//#define DSP_DEBUG

#include <Solver/DualDecomp/DdMWPara.h>

DdMWPara::DdMWPara(
		MPI_Comm          comm,   /**< MPI communicator */
		DdMaster *        master, /**< master problem */
		vector<DdWorker*> worker  /**< worker for finding lower bounds */):
DdMW(master, worker),
comm_(comm), sync_(false),
subcomm_(MPI_COMM_NULL), subcomm_size_(0),
num_comm_colors_(0), comm_color_(0), comm_key_(-1), comm_root_keys_(NULL), //comm_root_key_(-1),
lb_comm_(MPI_COMM_NULL), lb_group_(MPI_GROUP_EMPTY),
cg_comm_(MPI_COMM_NULL), cg_group_(MPI_GROUP_EMPTY), cg_comm_rank_(-1),
ub_comm_(MPI_COMM_NULL), ub_group_(MPI_GROUP_EMPTY), ub_comm_rank_(-1),
nsubprobs_(NULL), subprob_indices_(NULL), subprob_displs_(NULL)//, splitWorkers_(false)
{
	MPI_Comm_rank(comm, &comm_rank_);
	MPI_Comm_size(comm, &comm_size_);
}

DdMWPara::~DdMWPara() {
	FREE_ARRAY_PTR(comm_root_keys_);
	FREE_ARRAY_PTR(nsubprobs_);
	FREE_ARRAY_PTR(subprob_indices_);
	FREE_ARRAY_PTR(subprob_displs_);
}

DSP_RTN_CODE DdMWPara::init() {
	BGN_TRY_CATCH

	DdMW::init();

	int narrprocidx  = par_->getIntPtrParamSize("ARR_PROC_IDX"); /**< size of pointer parameter ARR_PROC_IDX */
	int * arrprocidx = par_->getIntPtrParam("ARR_PROC_IDX");     /**< pointer parameter ARR_PROC_IDX */

	if (comm_rank_ == 0)
	{
		nsubprobs_      = new int [comm_size_];
		subprob_displs_ = new int [comm_size_];
	}

	/** gather number of subproblems for each process */
	MPI_Gather(&narrprocidx, 1, MPI_INT, nsubprobs_, 1, MPI_INT, 0, comm_);

	if (comm_rank_ == 0)
	{
		int sum_of_nsubprobs = 0;
		for (int i = 0; i < comm_size_; ++i)
		{
			sum_of_nsubprobs += nsubprobs_[i];
			subprob_displs_[i] = i == 0 ? 0 : subprob_displs_[i-1] + nsubprobs_[i-1];
		}
		subprob_indices_ = new int [sum_of_nsubprobs];
	}

	/** gather subproblem indices for each process */
	MPI_Gatherv(arrprocidx, narrprocidx, MPI_INT, subprob_indices_, nsubprobs_, subprob_displs_, MPI_INT, 0, comm_);

#ifdef DSP_DEBUG
	if (comm_rank_ == 0)
		DSPdebugMessage("comm_size_ %d NumSubproblems %d\n", comm_size_, model_->getNumSubproblems());
#endif

	/** coloring the processors */
	if (comm_rank_ == 0)
		comm_color_ = comm_color_main;
	else
	{
		comm_color_ = comm_color_main;
		if (comm_rank_ > model_->getNumSubproblems() + 1)
			comm_color_ = comm_color_cg;
		if (comm_rank_ > 2 * model_->getNumSubproblems() + 1)
			comm_color_ = comm_color_ub;
	}

	/** number of rank colors */
	num_comm_colors_ = 1;
	if (comm_size_ > model_->getNumSubproblems() + 1)
		num_comm_colors_++;
	if (comm_size_ > 2 * model_->getNumSubproblems() + 1)
		num_comm_colors_++;

	/** split communication */
	MPI_Comm_split(comm_, comm_color_, comm_rank_, &subcomm_);
	MPI_Comm_size(subcomm_, &subcomm_size_);
	MPI_Comm_rank(subcomm_, &comm_key_);

	DSPdebugMessage("comm_rank_ %d comm_color_ %d subcomm_size_ %d comm_key_ %d\n",
			comm_rank_, comm_color_, subcomm_size_, comm_key_);

	/** allocate memory */
	comm_root_keys_ = new int [3];
	comm_root_keys_[0] = 0; /**< the root */
	comm_root_keys_[1] = 1;
	comm_root_keys_[2] = 1;
	/**< root rank for comm_color_cg */
	if (comm_size_ > model_->getNumSubproblems() + 1)
		comm_root_keys_[1] = model_->getNumSubproblems() + 1;
	/**< root rank for comm_color_ub */
	if (comm_size_ > 2 * model_->getNumSubproblems() + 1)
		comm_root_keys_[2] = 2 * model_->getNumSubproblems() + 1;

	/** MPI group */
	MPI_Group world_group;
	MPI_Comm_group(comm_, &world_group);
	int * ranks = NULL;

	/** create a lower bouding group */
	if (comm_size_ < model_->getNumSubproblems() + 1)
	{
		ranks = new int [comm_size_ - 1];
		for (int i = 1; i < comm_size_; ++i)
			ranks[i - 1] = i;
		MPI_Group_incl(world_group, comm_size_ - 1, ranks, &lb_group_);
		FREE_ARRAY_PTR(ranks);
	}
	else
	{
		ranks = new int [model_->getNumSubproblems()];
		for (int i = 1; i < model_->getNumSubproblems() + 1; ++i)
			ranks[i - 1] = i;
		MPI_Group_incl(world_group, model_->getNumSubproblems(), ranks, &lb_group_);
		FREE_ARRAY_PTR(ranks);
	}
	MPI_Comm_create_group(comm_, lb_group_, 0, &lb_comm_);

	/** create a cut generation group */
	if (parFeasCuts_ >= 0 || parOptCuts_ >= 0)
	{
		int nranks = 0;
		if (comm_size_ <= model_->getNumSubproblems() + 1)
		{
			nranks = comm_size_ - 1;
			ranks = new int [nranks];
			for (int i = 0; i < nranks; ++i)
				ranks[i] = i + 1;
		}
		else if (comm_size_ <= 2 * model_->getNumSubproblems() + 1)
		{
			nranks = comm_size_ - model_->getNumSubproblems() - 1;
			ranks = new int [nranks];
			for (int i = 0; i < nranks; ++i)
				ranks[i] = i + model_->getNumSubproblems() + 1;
		}
		else
		{
			nranks = model_->getNumSubproblems();
			ranks = new int [nranks];
			for (int i = 0; i < nranks; ++i)
				ranks[i] = i + model_->getNumSubproblems() + 1;
		}
		MPI_Group_incl(world_group, nranks, ranks, &cg_group_);
		FREE_ARRAY_PTR(ranks);

		MPI_Comm_create_group(comm_, cg_group_, 0, &cg_comm_);
		if (cg_comm_ != MPI_COMM_NULL)
		{
			MPI_Comm_rank(cg_comm_, &cg_comm_rank_);
			DSPdebugMessage("comm_ %d cg_group_ %d cg_comm_ %d cg_comm_rank_ %d\n", comm_, cg_group_, cg_comm_, cg_comm_rank_);
		}
	}

	/** create a upper bound group */
	if (parEvalUb_ >= 0)
	{
		int nranks = 0;
		if (comm_size_ <= model_->getNumSubproblems() + 1)
		{
			nranks = comm_size_ - 1;
			ranks = new int [nranks];
			for (int i = 0; i < nranks; ++i)
				ranks[i] = i + 1;
		}
		else if (comm_size_ <= 2 * model_->getNumSubproblems() + 1)
		{
			nranks = comm_size_ - model_->getNumSubproblems() - 1;
			ranks = new int [nranks];
			for (int i = 0; i < nranks; ++i)
				ranks[i] = i + model_->getNumSubproblems() + 1;
		}
		else if (comm_size_ <= 3 * model_->getNumSubproblems() + 1)
		{
			nranks = comm_size_ - 2 * model_->getNumSubproblems() - 1;
			ranks = new int [nranks];
			for (int i = 0; i < nranks; ++i)
				ranks[i] = i + 2 * model_->getNumSubproblems() + 1;
		}
		else
		{
			nranks = model_->getNumSubproblems();
			ranks = new int [nranks];
			for (int i = 0; i < nranks; ++i)
				ranks[i] = i + 2 * model_->getNumSubproblems() + 1;
		}
		MPI_Group_incl(world_group, nranks, ranks, &ub_group_);
		FREE_ARRAY_PTR(ranks);

		MPI_Comm_create_group(comm_, ub_group_, 0, &ub_comm_);
		if (ub_comm_ != MPI_COMM_NULL)
		{
			MPI_Comm_rank(ub_comm_, &ub_comm_rank_);
			DSPdebugMessage("comm_ %d ub_group_ %d ub_comm_ %d ub_comm_rank_ %d\n", comm_, ub_group_, ub_comm_, ub_comm_rank_);
		}
	}

	MPI_Group_free(&world_group);

#if 0
	if (comm_size_ > model_->getNumSubproblems() + 1)
	{
		splitWorkers_ = true;
		/** TODO could have many processors */
		if (comm_rank_ > model_->getNumSubproblems())
			comm_color_ = 1;
		MPI_Comm_split(comm_, comm_color_, comm_rank_, &subcomm_);
		MPI_Comm_size(subcomm_, &subcomm_size_);
		MPI_Comm_rank(subcomm_, &comm_key_);
		if (comm_key_ == 0)
		{
			MPI_Status status;
			if (comm_color_ == 0)
				MPI_Recv(&comm_root_key_, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, comm_, &status);
			else
				MPI_Send(&comm_rank_, 1, MPI_INT, 0, 0, comm_);
		}
#ifdef DSP_DEBUG
		for (int i = 0; i < comm_size_; ++i)
		{
			if (i == comm_rank_)
				DSPdebugMessage("comm_color_ %d subcomm_size_ %d comm_key %d\n", comm_color_, subcomm_size_, comm_key_);
			MPI_Barrier(comm_);
		}
#endif
	}
#endif

	/** free memory */
	arrprocidx = NULL;

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMWPara::finalize() {
#define FREE_MPI_COMM(_comm)      \
	if (_comm != MPI_COMM_NULL) { \
		MPI_Comm_free(&_comm);    \
		_comm = MPI_COMM_NULL;    \
	}
#define FREE_MPI_GROUP(_group)       \
	if (_group != MPI_GROUP_EMPTY) { \
		MPI_Group_free(&_group);     \
		_group = MPI_GROUP_EMPTY;    \
	}

	BGN_TRY_CATCH

	DdMW::finalize();

	FREE_MPI_GROUP(lb_group_)
	FREE_MPI_COMM(lb_comm_)
	FREE_MPI_GROUP(cg_group_)
	FREE_MPI_COMM(cg_comm_)
	FREE_MPI_GROUP(ub_group_)
	FREE_MPI_COMM(ub_comm_)
	FREE_MPI_COMM(subcomm_)

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}
