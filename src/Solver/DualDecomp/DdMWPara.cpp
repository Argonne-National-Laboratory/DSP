/*
 * DdMWPara.cpp
 *
 *  Created on: Apr 8, 2016
 *      Author: kibaekkim
 */

// #define DSP_DEBUG

#include "Model/TssModel.h"
#include "DdMWPara.h"

DdMWPara::DdMWPara(
		MPI_Comm     comm,   /**< MPI communicator */
		DecModel *   model,  /**< model pointer */
		DspParams *  par,    /**< parameters */
		DspMessage * message /**< message pointer */):
DdMW(model,par,message),
comm_(comm),
sync_(false),
subcomm_(MPI_COMM_NULL),
subcomm_group_(MPI_GROUP_NULL),
subcomm_size_(0),
subcomm_rank_(-1),
lb_comm_(MPI_COMM_NULL),
lb_group_(MPI_GROUP_NULL),
lb_comm_size_(0),
lb_comm_rank_(-1),
lb_comm_root_(-1),
has_cgub_comm_(false),
cgub_comm_(MPI_COMM_NULL),
cgub_group_(MPI_GROUP_NULL),
cgub_comm_size_(0),
cgub_comm_rank_(-1),
cgub_comm_root_(-1),
nsubprobs_(NULL),
subprob_indices_(NULL),
subprob_displs_(NULL) {
	MPI_Comm_rank(comm, &comm_rank_);
	MPI_Comm_size(comm, &comm_size_);
}

DdMWPara::DdMWPara(const DdMWPara& rhs) :
DdMW(rhs),
comm_(rhs.comm_),
comm_rank_(rhs.comm_rank_),
comm_size_(rhs.comm_size_),
sync_(rhs.sync_),
subcomm_(rhs.subcomm_),
subcomm_group_(rhs.subcomm_group_),
subcomm_size_(rhs.subcomm_size_),
subcomm_rank_(rhs.subcomm_rank_),
lb_comm_(rhs.lb_comm_),
lb_group_(rhs.lb_group_),
lb_comm_size_(rhs.lb_comm_size_),
lb_comm_rank_(rhs.lb_comm_rank_),
lb_comm_root_(rhs.lb_comm_root_),
has_cgub_comm_(rhs.has_cgub_comm_),
cgub_comm_(rhs.cgub_comm_),
cgub_group_(rhs.cgub_group_),
cgub_comm_size_(rhs.cgub_comm_size_),
cgub_comm_rank_(rhs.cgub_comm_rank_),
cgub_comm_root_(rhs.cgub_comm_root_) {
	// copy number of subproblems for each proc
	nsubprobs_ = new int [comm_size_];
	CoinCopyN(rhs.nsubprobs_, comm_size_, nsubprobs_);

	// copy subproblem indices for each proc
	int sum_of_nsubprobs = 0;
	for (int i = 0; i < comm_size_; ++i)
		sum_of_nsubprobs += nsubprobs_[i];
	subprob_indices_ = new int [sum_of_nsubprobs];
	CoinCopyN(rhs.subprob_indices_, sum_of_nsubprobs, subprob_indices_);

	// copy subproblem index locations for each proc
	subprob_displs_ = new int [comm_size_];
	CoinCopyN(rhs.subprob_displs_, comm_size_, subprob_displs_);
}

DdMWPara::~DdMWPara() {
	FREE_ARRAY_PTR(nsubprobs_);
	FREE_ARRAY_PTR(subprob_indices_);
	FREE_ARRAY_PTR(subprob_displs_);
}

DSP_RTN_CODE DdMWPara::init() {

	BGN_TRY_CATCH

	DdMW::init();

	DSPdebugMessage("Beginning of Processor Management.\n");

	/** create cgub_comm_ only for asynchronous computing */
	if (sync_ == false && (parFeasCuts_ >= 0 || parOptCuts_ >= 0 || parEvalUb_ >= 0))
		has_cgub_comm_ = true;

	/** TODO At least three processors are required for running CG and/or UB */
	if (has_cgub_comm_ == true && comm_size_ < 3)
	{
		if (parFeasCuts_ >= 0 || parOptCuts_ >= 0)
		{
			if (comm_rank_ == 0)
				message_->print(0, "Warning: Cut generation is disabled. Enabling cut generation requires at least three processors.\n");
			parFeasCuts_ = -1;
			parOptCuts_ = -1;
		}
		if (parEvalUb_ >= 0)
		{
			if (comm_rank_ == 0)
				message_->print(0, "Warning: Finding upper bounds is disabled. Enabling upper bounds requires at least three processors.\n");
			parEvalUb_ = -1;
		}
		has_cgub_comm_ = false;
	}
#ifdef DSP_DEBUG
	if (comm_rank_ == 0)
		DSPdebugMessage("comm_size_ %d NumSubproblems %d\n", comm_size_, model_->getNumSubproblems());
#endif

	/** create communication groups */
	createGroups();

	/** set subproblem indices to each processor */
	setSubproblemIndices();

	/** assign root key for each communication group */
	setRootKeys();

	DSPdebugMessage("End of Processor Management.\n");

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
	if (_group != MPI_GROUP_NULL) { \
		MPI_Group_free(&_group);     \
		_group = MPI_GROUP_NULL;    \
	}

	BGN_TRY_CATCH

	FREE_MPI_GROUP(lb_group_)
	FREE_MPI_COMM(lb_comm_)
	FREE_MPI_GROUP(cgub_group_)
	FREE_MPI_COMM(cgub_comm_)
	FREE_MPI_GROUP(subcomm_group_)
	FREE_MPI_COMM(subcomm_)

	DdMW::finalize();

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMWPara::createGroups() {

	MPI_Group world_group;
	int nranks;
	std::vector<int> ranks;

	BGN_TRY_CATCH

	/** get world group */
	MPI_Comm_group(comm_, &world_group);

	if (comm_size_ >= 2 * model_->getNumSubproblems() + 1)
	{
		/** Master-LB group */
		nranks = model_->getNumSubproblems() + 1;
		ranks.resize(nranks);
		for (int i = 0; i < nranks; ++i)
			ranks[i] = i;
		MPI_Group_incl(world_group, nranks, ranks.data(), &subcomm_group_);
		MPI_Comm_create_group(comm_, subcomm_group_, DSP_MPI_TAG_GROUP_SUB, &subcomm_);
		if (subcomm_ != MPI_COMM_NULL)
		{
			MPI_Comm_size(subcomm_, &subcomm_size_);
			MPI_Comm_rank(subcomm_, &subcomm_rank_);
		}
		/** LB group */
		MPI_Group_incl(world_group, nranks-1, ranks.data()+1, &lb_group_);
		MPI_Comm_create_group(comm_, lb_group_, DSP_MPI_TAG_GROUP_LB, &lb_comm_);
		if (lb_comm_ != MPI_COMM_NULL)
		{
			MPI_Comm_rank(lb_comm_, &lb_comm_rank_);
			MPI_Comm_size(lb_comm_, &lb_comm_size_);
		}
		/** CG-UB group */
		if (has_cgub_comm_)
		{
			nranks = model_->getNumSubproblems();
			ranks.resize(nranks);
			for (int i = 0; i < nranks; ++i)
				ranks[i] = i+model_->getNumSubproblems()+1;
			MPI_Group_incl(world_group, nranks, ranks.data(), &cgub_group_);
			MPI_Comm_create_group(comm_, cgub_group_, DSP_MPI_TAG_GROUP_CGUB, &cgub_comm_);
			if (cgub_comm_ != MPI_COMM_NULL)
			{
				MPI_Comm_rank(cgub_comm_, &cgub_comm_rank_);
				MPI_Comm_size(cgub_comm_, &cgub_comm_size_);
			}
		}
	}
	else
	{
		/** TODO equally split the processors for now */
		if (has_cgub_comm_) {
			//nranks = ceil(1.0*(comm_size_-1)/2) + 1;
			nranks = comm_size_ * par_->getDblParam("DD/WORKER_RATIO");
		} else
			nranks = min(comm_size_,model_->getNumSubproblems() + 1);
		ranks.resize(nranks);
		for (int i = 0; i < nranks; ++i)
			ranks[i] = i;
		int ierr = MPI_Group_incl(world_group, nranks, ranks.data(), &subcomm_group_);
		ierr = MPI_Comm_create_group(comm_, subcomm_group_, DSP_MPI_TAG_GROUP_SUB, &subcomm_);
		if (subcomm_ != MPI_COMM_NULL)
		{
			MPI_Comm_size(subcomm_, &subcomm_size_);
			MPI_Comm_rank(subcomm_, &subcomm_rank_);
			DSPdebugMessage("Rank %d: subcomm_size_ %d, subcomm_rank_ %d\n",
					comm_rank_, subcomm_size_, subcomm_rank_);
		}
		/** LB group */
		DSPdebug(DspMessage::printArray(nranks-1, ranks.data()+1));
		MPI_Group_incl(world_group, nranks-1, ranks.data()+1, &lb_group_);
		MPI_Comm_create_group(comm_, lb_group_, DSP_MPI_TAG_LB, &lb_comm_);
		if (lb_comm_ != MPI_COMM_NULL)
		{
			MPI_Comm_rank(lb_comm_, &lb_comm_rank_);
			MPI_Comm_size(lb_comm_, &lb_comm_size_);
		}
		/** CG-UB group */
		if (has_cgub_comm_)
		{
			MPI_Group_excl(world_group, nranks, ranks.data(), &cgub_group_);
			MPI_Comm_create_group(comm_, cgub_group_, DSP_MPI_TAG_CGUB, &cgub_comm_);
			if (cgub_comm_ != MPI_COMM_NULL)
			{
				MPI_Comm_rank(cgub_comm_, &cgub_comm_rank_);
				MPI_Comm_size(cgub_comm_, &cgub_comm_size_);
			}
		}
	}
#ifdef DSP_DEBUG
	for (int i = 0; i < comm_size_; ++i)
	{
		if (i == comm_rank_)
			DSPdebugMessage("comm_rank_ %d, subcomm_rank_ %d, lb_comm_rank_ %d, cgub_comm_rank_ %d, cgub_comm_ %s\n",
					comm_rank_, subcomm_rank_, lb_comm_rank_, cgub_comm_rank_, cgub_comm_ == MPI_COMM_NULL ? "null" : "created");
		MPI_Barrier(comm_);
	}
#endif

	/** free world group */
	MPI_Group_free(&world_group);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

/** TODO This may not be the best place for this function. */
DSP_RTN_CODE DdMWPara::generateBendersCuts(
		MPI_Comm  comm,      /**< communicator */
		int       comm_rank, /**< rank in comm */
		Solutions solutions, /**< solutions at which cuts are generated */
		OsiCuts & cuts       /**< cuts generated */) {
#define FREE_MEMORY \
	FREE_ARRAY_PTR(aggcut)

	int ret = DSP_STAT_MW_CONTINUE;

#ifdef DSP_HAS_SCIP
	if (solutions.size() == 0) return ret;
	if (parFeasCuts_ < 0 && parOptCuts_ < 0) return ret;
	if (model_->isStochastic() == false)
	{
		message_->print(0, "This problem is not a stochastic program. Benders cut option is valid only for stochastic programming.\n");
		parFeasCuts_ = -1;
		parOptCuts_ = -1;
		return ret;
	}

	vector<int> cuttype;
	TssModel * tssmodel = NULL;
	double * aggcut = NULL;
	double aggrhs;

	BGN_TRY_CATCH

	/** retrieve DdWorkerCGBd */
	DdWorkerCGBd * workercg = NULL;
	for (unsigned i = 0; i < worker_.size(); ++i)
	{
		if (worker_[i]->getType() == DdWorker::CGBd)
		{
			workercg = dynamic_cast<DdWorkerCGBd*>(worker_[i]);
			DSPdebugMessage("Rank %d works for generating Benders cuts.\n", comm_rank_);
			break;
		}
	}
	if (workercg == NULL)
	{
		DSPdebugMessage("Rank %d does not work for generating Benders cuts.\n", comm_rank_);
		return ret;
	}

	/** timing */
	double cputime = CoinCpuTime();
	double walltime = CoinGetTimeOfDay();

	/** downcast to TssModel */
	tssmodel = dynamic_cast<TssModel*>(model_);

	/** allocate memory */
	aggcut = new double [tssmodel->getNumCols(0) + 1];

	/** generate cuts */
	for (unsigned i = 0; i < solutions.size(); ++i)
	{
		/** create local cut pools */
		OsiCuts localcuts, allcuts;

		/** generate cuts at a solution */
		DSP_RTN_CHECK_THROW(workercg->generateCuts(solutions[i], &localcuts, cuttype));
		DSPdebugMessage("Rank %d: Benders cut generator generated %d cuts for solution %d.\n",
				comm_rank_, localcuts.sizeCuts(), i);

		/** check if there exists a feasibility cut */
		bool hasFeasibilityCuts = false;
		bool hasOptimalityCuts = true;
		for (unsigned j = 0; j < cuttype.size(); ++j)
		{
			if (cuttype[j] == DdWorkerCGBd::Feas)
			{
				hasFeasibilityCuts = true;
				hasOptimalityCuts = false;
				break;
			}
			if (cuttype[j] != DdWorkerCGBd::Opt)
				hasOptimalityCuts = false;
		}

		/** MPI_Allreduce */
		bool genFeasibilityCuts = false;
		bool genOptimalityCuts  = false;
		MPI_Allreduce(&hasFeasibilityCuts, &genFeasibilityCuts, 1, MPI_C_BOOL, MPI_LOR, comm);
		MPI_Allreduce(&hasOptimalityCuts, &genOptimalityCuts, 1, MPI_C_BOOL, MPI_LAND, comm);
		DSPdebugMessage2("Rank %d: hasFeasibilityCuts = %s, hasOptimalityCuts = %s\n",
				comm_rank_, hasFeasibilityCuts ? "true" : "false", hasOptimalityCuts ? "true" : "false");

		if ((genFeasibilityCuts == true && parFeasCuts_ >= 0) ||
			(genOptimalityCuts == true && parOptCuts_ >= 0))
		{
			MPIgatherOsiCuts(comm, localcuts, allcuts);
			/** free local cuts */
			for (int j = 0; j < localcuts.sizeCuts(); ++j)
			{
				OsiRowCut * rc = localcuts.rowCutPtr(j);
				FREE_PTR(rc);
			}
			localcuts.dumpCuts();

			if (comm_rank == 0)
			{
				DSPdebugMessage("Rank %d has %d cuts gathered.\n", comm_rank_, allcuts.sizeCuts());
				if (genOptimalityCuts)
				{
					/** construct optimality cut */
					CoinZeroN(aggcut, tssmodel->getNumCols(0) + 1);
					for (int j = 0; j < tssmodel->getNumCols(0); ++j)
						aggcut[j] = -(tssmodel->getObjCore(0)[j]);
					aggrhs = 0.0;

					for (int j = 0; j < allcuts.sizeCuts(); ++j)
					{
						OsiRowCut * rc = allcuts.rowCutPtr(j);
						CoinPackedVector cutrow = rc->row();
						for (int k = 0; k < cutrow.getNumElements(); ++k)
							aggcut[cutrow.getIndices()[k]] += cutrow.getElements()[k];
						aggrhs += rc->lb();
						DSPdebugMessage2("allcuts[%d]:\n", j);
						DSPdebug2(rc->print());
					}
					DSPdebug2(message_->printArray(tssmodel->getNumCols(0) + 1, aggcut));

					CoinPackedVector orow;
					for (int j = 0; j < tssmodel->getNumCols(0); ++j)
						if (fabs(aggcut[j]) > 1.0e-8)
							orow.insert(j, aggcut[j]);
					orow.insert(tssmodel->getNumCols(0) + tssmodel->getNumCols(1), 1.0);

					OsiRowCut ocut;
					ocut.setRow(orow);
					ocut.setUb(COIN_DBL_MAX);
					ocut.setLb(aggrhs);

					DSPdebug2(ocut.print());
					cuts.insertIfNotDuplicate(ocut);
				}

				if (genFeasibilityCuts)
				{
					/** copy cut pointers */
					for (int j = 0; j < allcuts.sizeCuts(); ++j)
					{
						OsiRowCut * rc = allcuts.rowCutPtr(j);
						DSPdebug(rc->print());
						cuts.insertIfNotDuplicate(*rc);
						ret = DSP_STAT_MW_RESOLVE;
					}
				}
			}
		}

		/** free cuts */
		for (int j = 0; j < allcuts.sizeCuts(); ++j)
		{
			OsiRowCut * rc = allcuts.rowCutPtr(j);
			FREE_PTR(rc);
		}
		allcuts.dumpCuts();
	}

	/** log timing results */
	workercg->s_statuses_.push_back(DSP_RTN_OK);
	workercg->s_primobjs_.push_back(0.0);
	workercg->s_dualobjs_.push_back(0.0);
	workercg->s_cputimes_.push_back(CoinCpuTime() - cputime);
	workercg->s_walltimes_.push_back(CoinGetTimeOfDay() - walltime);

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY
#endif

	return ret;
#undef FREE_MEMORY
}

DSP_RTN_CODE DdMWPara::setSubproblemIndices() {
#define FREE_MEMORY \
	FREE_ARRAY_PTR(arrprocidx)

	int narrprocidx = 0;
	int * arrprocidx = NULL;

	BGN_TRY_CATCH

	/** allocate memory */
	nsubprobs_      = new int [comm_size_];
	subprob_displs_ = new int [comm_size_];

	/** assign subproblems to each processor */
	if (lb_comm_ != MPI_COMM_NULL)
	{
		narrprocidx = 0;
		for (int i = lb_comm_rank_; i < model_->getNumSubproblems(); i += lb_comm_size_)
			narrprocidx++;
		arrprocidx = new int [narrprocidx];
		for (int i = lb_comm_rank_, j = 0; i < model_->getNumSubproblems(); i += lb_comm_size_)
			arrprocidx[j++] = i;
	}
	else if (cgub_comm_ != MPI_COMM_NULL)
	{
		narrprocidx = 0;
		for (int i = cgub_comm_rank_; i < model_->getNumSubproblems(); i += cgub_comm_size_)
			narrprocidx++;
		arrprocidx = new int [narrprocidx];
		for (int i = cgub_comm_rank_, j = 0; i < model_->getNumSubproblems(); i += cgub_comm_size_)
			arrprocidx[j++] = i;
	}

#ifdef DSP_DEBUG
	for (int i = 0; i < comm_size_; ++i)
	{
		if (i == comm_rank_)
		{
			DSPdebugMessage("Rank %d has the following subproblem indices (%d):\n", comm_rank_, narrprocidx);
			DspMessage::printArray(narrprocidx, arrprocidx);
		}
		MPI_Barrier(comm_);
	}
#endif

	/** gather number of subproblems for each process */
	MPI_Allgather(&narrprocidx, 1, MPI_INT, nsubprobs_, 1, MPI_INT, comm_);

	/** gather subproblem indices for each process */
	int sum_of_nsubprobs = 0;
	for (int i = 0; i < comm_size_; ++i)
	{
		sum_of_nsubprobs += nsubprobs_[i];
		subprob_displs_[i] = i == 0 ? 0 : subprob_displs_[i-1] + nsubprobs_[i-1];
	}
	subprob_indices_ = new int [sum_of_nsubprobs];
	MPI_Allgatherv(arrprocidx, narrprocidx, MPI_INT, subprob_indices_, nsubprobs_, subprob_displs_, MPI_INT, comm_);
#ifdef DSP_DEBUG
	if (comm_rank_ == 0)
	{
		DSPdebugMessage("Rank %d has subprob_indices_:\n", comm_rank_);
		DspMessage::printArray(sum_of_nsubprobs, subprob_indices_);
	}
#endif

	/** set parameters */
	par_->setIntPtrParamSize("ARR_PROC_IDX", narrprocidx);
	for (int i = 0; i < narrprocidx; ++i)
		par_->setIntPtrParam("ARR_PROC_IDX", i, arrprocidx[i]);

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

DSP_RTN_CODE DdMWPara::setRootKeys() {

	BGN_TRY_CATCH

	MPI_Status status;

	if (comm_rank_ == 0)
	{
		MPI_Recv(&lb_comm_root_, 1, MPI_INT, MPI_ANY_SOURCE, DSP_MPI_TAG_LB, comm_, &status);
		if (has_cgub_comm_)
			MPI_Recv(&cgub_comm_root_, 1, MPI_INT, MPI_ANY_SOURCE, DSP_MPI_TAG_CGUB, comm_, &status);
	}
	if (lb_comm_rank_ == 0)
		MPI_Send(&comm_rank_, 1, MPI_INT, 0, DSP_MPI_TAG_LB, comm_);
	if (cgub_comm_rank_ == 0)
		MPI_Send(&comm_rank_, 1, MPI_INT, 0, DSP_MPI_TAG_CGUB, comm_);
	MPI_Bcast(&lb_comm_root_, 1, MPI_INT, 0, comm_);
	MPI_Bcast(&lb_comm_size_, 1, MPI_INT, lb_comm_root_, comm_);
	if (has_cgub_comm_)
	{
		MPI_Bcast(&cgub_comm_root_, 1, MPI_INT, 0, comm_);
		MPI_Bcast(&cgub_comm_size_, 1, MPI_INT, cgub_comm_root_, comm_);
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}
