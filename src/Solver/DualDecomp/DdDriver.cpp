/*
 * DdDriver.cpp
 *
 *  Created on: Feb 16, 2016
 *      Author: kibaekkim
 */

#include "Solver/DualDecomp/DdDriver.h"
#include "Solver/DualDecomp/DdMasterTr.h"
#include "Solver/DualDecomp/DdMasterDsb.h"
#include "Solver/DualDecomp/DdMasterSubgrad.h"
#include "Solver/DualDecomp/DdMasterReg.h"

DdDriver::DdDriver(DspParams * par, DecModel * model):
	DspDriver(par, model),
	comm_(-1), comm_rank_(0), comm_size_(0),
	mw_(NULL), master_(NULL), worker_(NULL),
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

	/** create worker */
	worker_ = new DdWorker(par_, model_, message_);

	END_TRY_CATCH(;)
}


/** constructor with MPI */
DdDriver::DdDriver(DspParams * par, DecModel * model, MPI_Comm comm):
	DspDriver(par, model), comm_(comm),
	mw_(NULL), master_(NULL), worker_(NULL),
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
		worker_ = new DdWorker(par_, model_, message_);
	}

	END_TRY_CATCH(;)
}

DdDriver::~DdDriver()
{
	FREE_PTR(mw_);
	FREE_PTR(master_);
	FREE_PTR(worker_);
}

/** initialize */
STO_RTN_CODE DdDriver::init()
{
	BGN_TRY_CATCH

	/** initialize master/worker */
	if (master_) master_->init();
	if (worker_) worker_->init();

	/** create master/worker framework */
	mw_ = new DdMW(comm_, master_, worker_);

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
}

/** run */
STO_RTN_CODE DdDriver::run()
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

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
}
