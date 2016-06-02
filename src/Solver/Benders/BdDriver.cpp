/*
 * BdDriver.cpp
 *
 *  Created on: Feb 17, 2016
 *      Author: kibaekkim
 */

#include "Solver/Benders/BdDriver.h"
#include "Solver/Benders/SCIPconshdlrBenders.h"
#include "Solver/Benders/SCIPconshdlrBendersWorker.h"
#include "Solver/DualDecomp/DdDriver.h"
#include "SolverInterface/SolverInterfaceScip.h"

BdDriver::BdDriver(DspParams * par, DecModel * model):
	DspDriver(par, model),
	useMPI_(false), comm_(MPI_UNDEFINED), comm_rank_(0), comm_size_(0),
	mw_(NULL), master_(NULL), worker_(NULL),
	numPriorities_(0), priorities_(NULL)
{
	/** create master and worker */
	master_ = new BdMaster(par_, model_, message_);
	worker_ = new BdWorker(par_, model_, message_);
}

BdDriver::BdDriver(DspParams * par, DecModel * model, MPI_Comm comm):
	DspDriver(par, model),
	useMPI_(true), comm_(comm),
	mw_(NULL), master_(NULL), worker_(NULL),
	numPriorities_(0), priorities_(NULL)
{
	MPI_Comm_rank(comm_, &comm_rank_);
	MPI_Comm_size(comm_, &comm_size_);

	/** create master and worker */
	if (comm_rank_ == 0)
		master_ = new BdMaster(par_, model_, message_);
	else
		worker_ = new BdWorker(par_, model_, message_);
}

BdDriver::~BdDriver()
{
	FREE_PTR(mw_);
	FREE_PTR(master_);
	FREE_PTR(worker_);
}

/** initialize */
DSP_RTN_CODE BdDriver::init()
{
	BGN_TRY_CATCH

	/** initialize master */
	if (master_)
	{
		master_->init();
		primsol_ = new double [model_->getFullModelNumCols()];
	}

	/** initialize worker */
	if (worker_) worker_->init();

	/** create master/worker framework */
	if (useMPI_)
		mw_ = new BdMW(comm_, master_, worker_);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

/** run */
DSP_RTN_CODE BdDriver::run()
{
	BGN_TRY_CATCH

	if (dualobj_ <= -COIN_DBL_MAX)
	{
		/** time stamp */
		double swtime = CoinGetTimeOfDay();

		/** find lower bound first if it is not given by user */
		message_->print(1, "Finding global lower bound ...");
		findLowerBound();
		message_->print(1, " (time elapsed: %.2f sec) -> Lower bound %e\n", CoinGetTimeOfDay() - swtime, dualobj_);
	}

	if (master_)
	{
		/** set constraint handler */
		master_->setConshdlr(constraintHandler());

		/** set dual bound, adding lower bounding constraint */
		master_->setDualObjective(dualobj_);

		/** set branching priorities */
		if (numPriorities_ > 0)
			master_->setBranchingPriority(numPriorities_, priorities_);

		/** set initial solutions */
		if (initsols_.size() > 0)
			master_->setSolutions(initsols_);
	}

	/** tic */
	cputime_  = CoinCpuTime();
	walltime_ = CoinGetTimeOfDay();

	/** run */
	if (!useMPI_)
	{
		/** solve */
		master_->solve();

		/** master solution */
		CoinCopyN(master_->getPrimalSolution(), model_->getNumSubproblemCouplingCols(0), primsol_);

		/** get subproblem objective values and solutions */
		int pos = model_->getNumSubproblemCouplingCols(0);
		for (int i = 0; i < model_->getNumSubproblems(); ++i)
		{
			int s = worker_->getBdSubPtr()->getSubprobIndex(i);
			CoinCopyN(worker_->getBdSubPtr()->getSolution(i),
					worker_->getBdSubPtr()->getNumCols(i),
					primsol_ + pos);
			pos += worker_->getBdSubPtr()->getNumCols(i);
		}
	}
	else
	{
		mw_->run();

		/** solution */
		CoinCopyN(mw_->getPrimalSolution(), model_->getFullModelNumCols(), primsol_);
	}

	/** toc */
	cputime_  = CoinCpuTime() - cputime_;
	walltime_ = CoinGetTimeOfDay() - walltime_;

	if (master_)
	{
		/** solution status */
		status_  = master_->getStatus();
		/** objective values */
		primobj_ = master_->getPrimalObjective();
		dualobj_ = master_->getDualObjective();
		/** number of nodes */
		numNodes_ = master_->getSiPtr()->getNumNodes();
		/** number of iterations */
		numIterations_ = master_->getSiPtr()->getIterationCount();
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE BdDriver::setAuxVarData(int size, double* obj, double* clbd, double* cubd)
{
	if (size <= 0) return DSP_RTN_ERR;
	if (!master_) return DSP_RTN_OK;

	BGN_TRY_CATCH

	master_->setAuxVarData(size, obj, clbd, cubd);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE BdDriver::setPriorities(
		int   size,      /**< size of array */
		int * priorities /**< branch priority */)
{
	if (size <= 0) return DSP_RTN_ERR;
	if (!master_) return DSP_RTN_OK;

	BGN_TRY_CATCH

	numPriorities_ = size;
	if (priorities_ == NULL)
		priorities_ = new int [numPriorities_];
	CoinCopyN(priorities, numPriorities_, priorities_);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE BdDriver::setSolution(
		int      size,    /**< size of array */
		double * solution /**< solution */)
{
	if (size <= 0) return DSP_RTN_ERR;
	if (!master_) return DSP_RTN_OK;

	BGN_TRY_CATCH

	initsols_.push_back(new CoinPackedVector(size, solution));

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

/** find lower bound */
DSP_RTN_CODE BdDriver::findLowerBound()
{
#define FREE_MEMORY \
	FREE_PTR(dd);

	DdDriver * dd = NULL;

	BGN_TRY_CATCH

	/** TODO use dual decomposition */
	dd = new DdDriver(par_, model_, comm_);
	dd->init();
	dd->run();

	/** set lower bound */
	dualobj_ = dd->getPrimalObjectiveValue();

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

SCIPconshdlrBenders * BdDriver::constraintHandler()
{
	SCIPconshdlrBenders * conshdlr = NULL;

	BGN_TRY_CATCH

	/** get solver interface */
	SolverInterfaceScip * si = dynamic_cast<SolverInterfaceScip*>(master_->getSiPtr());

	if (!useMPI_)
	{
		/** serial Benders */
		conshdlr = new SCIPconshdlrBenders(si->getSCIP(), par_->getIntParam("BD/CUT_PRIORITY"));
	}
	else
	{
		/** MPI Benders */
		conshdlr = new SCIPconshdlrBendersWorker(si->getSCIP(), par_->getIntParam("BD/CUT_PRIORITY"), comm_);
	}
	conshdlr->setNumSubprobs(model_->getNumSubproblems());
	conshdlr->setBdSub(worker_->getBdSubPtr());
	conshdlr->setOriginalVariables(si->getNumCols(), si->getSCIPvars(), par_->getIntParam("BD/NUM_CUTS_PER_ITER"));

	END_TRY_CATCH_RTN(;,NULL)

	return conshdlr;
}

DSP_RTN_CODE BdDriver::finalize() {
	return DSP_RTN_OK;
}
