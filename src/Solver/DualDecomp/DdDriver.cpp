/*
 * DdDriver.cpp
 *
 *  Created on: Feb 16, 2016
 *      Author: kibaekkim
 */

#define DSP_DEBUG

#include "Solver/DualDecomp/DdDriver.h"

/** constructor */
DdDriver::DdDriver(DspParams * par, DecModel * model):
DspDriver(par, model),
comm_(MPI_COMM_NULL),
comm_rank_(-1),
comm_size_(0),
mw_(NULL)
{}

/** constructor with MPI */
DdDriver::DdDriver(DspParams * par, DecModel * model, MPI_Comm comm):
DspDriver(par, model),
comm_(comm),
mw_(NULL)
{
	BGN_TRY_CATCH

	MPI_Comm_rank(comm_, &comm_rank_); /**< get process ID */
	MPI_Comm_size(comm_, &comm_size_); /**< get number of processes */

	END_TRY_CATCH(;)
}

DdDriver::~DdDriver()
{
	FREE_PTR(mw_);
}

/** initialize */
DSP_RTN_CODE DdDriver::init()
{
	BGN_TRY_CATCH
#if 1
	/** create Master-Worker framework */
	if (par_->getBoolParam("DD/ASYNC"))
		mw_ = new DdMWAsync(comm_, model_, par_, message_);
	else
		mw_ = new DdMWSync(comm_, model_, par_, message_);
	/** initialize master-worker framework */
	mw_->init();
#else
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

	/** create master/worker framework */
	if (par_->getBoolParam("DD/ASYNC"))
		mw_ = new DdMWAsync(comm_, master_, worker_);
	else
		mw_ = new DdMWSync(comm_, master_, worker_);
#endif

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

	/** retrieve master pointer */
	DdMaster * master = mw_->master_;
	if (master)
	{
		status_ = master->getStatus();
		primobj_ = master->getBestPrimalObjective();
		dualobj_ = master->getBestDualObjective();
		primsol_ = new double [master->getSiPtr()->getNumCols()];
		CoinCopyN(master->getPrimalSolution(), master->getSiPtr()->getNumCols(), primsol_);
		numNodes_ = master->getSiPtr()->getNumNodes();
		numIterations_ = master->getSiPtr()->getIterationCount();
	}
	/** nullify master pointer */
	master = NULL;

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdDriver::finalize()
{
	BGN_TRY_CATCH

	/** finalize master-worker framework */
	mw_->finalize();

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}
