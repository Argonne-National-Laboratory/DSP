/*
 * DdDriver.cpp
 *
 *  Created on: Feb 16, 2016
 *      Author: kibaekkim
 */

#define DSP_DEBUG

#include "Solver/DualDecomp/DdDriver.h"
#include "Solver/DualDecomp/DdMasterAtr.h"
#include "Solver/DualDecomp/DdMasterTr.h"
#include "Solver/DualDecomp/DdMasterDsb.h"
#include "Solver/DualDecomp/DdMasterSubgrad.h"
#include "Solver/DualDecomp/DdMasterReg.h"

DdDriver::DdDriver(DspParams * par, DecModel * model):
	DspDriver(par, model),
	comm_(-1), comm_rank_(0), comm_size_(0),
	mw_(NULL), master_(NULL),
	num_infeasible_solutions_(0)
{
	BGN_TRY_CATCH

	/**< number of workers */
	int nworkers = 0;
	/**< maximum number of subproblems considered in a worker */
	int maxnumsubprobs = model_->getNumSubproblems();

	/** create master */
	switch (par_->getIntParam("DD/MASTER_ALGO"))
	{
	case Simplex:
	case IPM:
	case IPM_Feasible:
		master_ = new DdMasterTr(par_, model_, message_, nworkers, maxnumsubprobs);
		break;
	case DSBM:
		master_ = new DdMasterDsb(par_, model_, message_, nworkers, maxnumsubprobs);
		break;
	case Subgradient:
		master_ = new DdMasterSubgrad(par_, model_, message_, nworkers, maxnumsubprobs);
		break;
	case Regularize_Bundle:
		master_ = new DdMasterReg(par_, model_, message_, nworkers, maxnumsubprobs);
		break;
	}

	/** create a lower bounder */
	worker_.push_back(new DdWorkerLB(par_, model_, message_));

	/** create a Benders cut generator */
	worker_.push_back(new DdWorkerCGBd(par_, model_, message_));

	/** create a upper bounder */
	worker_.push_back(new DdWorkerUB(par_, model_, message_));

	END_TRY_CATCH(;)
}


/** constructor with MPI */
DdDriver::DdDriver(DspParams * par, DecModel * model, MPI_Comm comm):
	DspDriver(par, model), comm_(comm),
	mw_(NULL), master_(NULL),
	num_infeasible_solutions_(0)
{
	BGN_TRY_CATCH

	MPI_Comm_rank(comm_, &comm_rank_); /**< get process ID */
	MPI_Comm_size(comm_, &comm_size_); /**< get number of processes */
	/**< maximum number of subproblems considered in a worker */
	int maxnumsubprobs = model_->getNumSubproblems();

	if (comm_rank_ == 0)
	{
		/** create master */
		switch (par_->getIntParam("DD/MASTER_ALGO"))
		{
		case Simplex:
		case IPM:
		case IPM_Feasible:
			if (par_->getBoolParam("DD/ASYNC"))
				master_ = new DdMasterAtr(par_, model_, message_, comm_size_, maxnumsubprobs);
			else
				master_ = new DdMasterTr(par_, model_, message_, comm_size_, maxnumsubprobs);
			break;
		case DSBM:
			master_ = new DdMasterDsb(par_, model_, message_, comm_size_, maxnumsubprobs);
			break;
		case Subgradient:
			master_ = new DdMasterSubgrad(par_, model_, message_, comm_size_, maxnumsubprobs);
			break;
		case Regularize_Bundle:
			master_ = new DdMasterReg(par_, model_, message_, comm_size_, maxnumsubprobs);
			break;
		}
	}
	else
	{
		/** create worker */
		if (comm_rank_ <= model_->getNumSubproblems())
		{
			DSPdebugMessage("Rank %d creates a worker for lower bounds.\n", comm_rank_);
			worker_.push_back(new DdWorkerLB(par_, model_, message_));
		}

		/** create worker for upper bounds */
		if (comm_rank_ > model_->getNumSubproblems() || comm_size_ - 1 <= model_->getNumSubproblems())
		{
			DSPdebugMessage("Rank %d creates a worker for Benders cut generation.\n", comm_rank_);
			worker_.push_back(new DdWorkerCGBd(par_, model_, message_));

			DSPdebugMessage("Rank %d creates a worker for upper bounds.\n", comm_rank_);
			worker_.push_back(new DdWorkerUB(par_, model_, message_));
		}
	}

	END_TRY_CATCH(;)
}

DdDriver::~DdDriver()
{
	FREE_PTR(mw_);
	FREE_PTR(master_);
	DSPdebugMessage("rank %d: worker size %lu\n", comm_rank_, worker_.size());
	for (unsigned i = 0; i < worker_.size(); ++i)
	{
		DSPdebugMessage("rank %d: free memory for worker type %d\n", comm_rank_, worker_[i]->getType());
		FREE_PTR(worker_[i]);
	}
}

/** initialize */
DSP_RTN_CODE DdDriver::init()
{
	BGN_TRY_CATCH

	/** initialize master/worker */
	if (master_) master_->init();
	for (unsigned i = 0; i < worker_.size(); ++i)
		worker_[i]->init();

	/** create master/worker framework */
	if (par_->getBoolParam("DD/ASYNC"))
		mw_ = new DdMWAsync(comm_, master_, worker_);
	else
		mw_ = new DdMWSync(comm_, master_, worker_);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

/** run */
DSP_RTN_CODE DdDriver::run()
{
	BGN_TRY_CATCH

	/** tic */
	cputime_  = CoinCpuTime();
	walltime_ = CoinGetTimeOfDay();

	/** run */
	mw_->run();

	/** toc */
	cputime_  = CoinCpuTime() - cputime_;
	walltime_ = CoinGetTimeOfDay() - walltime_;

	if (master_)
	{
		/** solution status */
		status_ = master_->getStatus();
		/** objective value */
		primobj_ = master_->getBestPrimalObjective();
		dualobj_ = master_->getBestDualObjective();
		/** solution */
		primsol_ = new double [master_->getSiPtr()->getNumCols()];
		CoinCopyN(master_->getPrimalSolution(), master_->getSiPtr()->getNumCols(), primsol_);
		/** number of nodes */
		numNodes_ = master_->getSiPtr()->getNumNodes();
		/** number of iterations */
		numIterations_ = master_->getSiPtr()->getIterationCount();
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdDriver::finalize()
{
	BGN_TRY_CATCH

	mw_->finalize();

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}
