/*
 * BdDriverMpi.cpp
 *
 *  Created on: Jul 10, 2016
 *      Author: kibaekkim
 */

//#define DSP_DEBUG
#include "Solver/Benders/BdDriverMpi.h"
#include "Solver/Benders/BdMWMpi.h"
#include "Solver/DualDecomp/DdDriverMpi.h"

BdDriverMpi::BdDriverMpi(
		DspParams* par,
		DecModel* model,
		MPI_Comm comm):
BdDriver(par, model),
comm_(comm)
{
	MPI_Comm_rank(comm, &comm_rank_);
	MPI_Comm_size(comm, &comm_size_);
}

BdDriverMpi::~BdDriverMpi() {
	// TODO Auto-generated destructor stub
}

DSP_RTN_CODE BdDriverMpi::init()
{
	BGN_TRY_CATCH

	if (comm_rank_ == 0)
		primsol_ = new double [model_->getFullModelNumCols()];

	/** create and initialize master-worker */
	mw_ = new BdMWMpi(comm_, model_, par_, message_);
	DSP_RTN_CHECK_THROW(mw_->init());

	if (comm_rank_ > 0)
		message_->logLevel_ = -999;

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE BdDriverMpi::run()
{
	BGN_TRY_CATCH

	/** tic */
	cputime_  = CoinCpuTime();
	walltime_ = CoinGetTimeOfDay();

	/** find lower bound first if it is not given by user */
	if (dualobj_ <= -COIN_DBL_MAX)
		DSP_RTN_CHECK_THROW(findLowerBound());

	if (comm_rank_ == 0)
	{
		/** set branching priorities */
		if (numPriorities_ > 0)
			DSP_RTN_CHECK_THROW(mw_->getMasterPtr()->setBranchingPriority(numPriorities_, priorities_));

		/** set initial solutions */
		if (initsols_.size() > 0)
			DSP_RTN_CHECK_THROW(mw_->getMasterPtr()->setSolutions(initsols_));
	}

	/** run */
	mw_->getMasterPtr()->setTimeLimit(
			par_->getDblParam("BD/WALL_LIM") - CoinGetTimeOfDay() + walltime_);
	DSP_RTN_CHECK_THROW(mw_->run());
	for (int i = 0; i < comm_size_; ++i)
	{
		if (i == comm_rank_)
			printf("Rank %d did run the Benders decomposition.\n", i);
		MPI_Barrier(comm_);
	}

	/** toc */
	cputime_  = CoinCpuTime() - cputime_;
	walltime_ = CoinGetTimeOfDay() - walltime_;

	/** collect solutions from the master and the workers */
	DSPdebugMessage("Rank %d: Collecting solutions...\n", comm_rank_);
	DSP_RTN_CHECK_THROW(collectSolution());

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE BdDriverMpi::finalize()
{
	BGN_TRY_CATCH

	mw_->finalize();

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE BdDriverMpi::findLowerBound()
{
#define FREE_MEMORY \
	FREE_PTR(dd);

	DdDriver * dd = NULL;

	BGN_TRY_CATCH

	message_->print(1, "Finding a good lower bound using Dual Decomposition...\n");

	/** TODO use dual decomposition */
	dd = new DdDriverMpi(par_, model_, comm_);
	DSP_RTN_CHECK_THROW(dd->init());
	DSP_RTN_CHECK_THROW(dd->run());
	DSP_RTN_CHECK_THROW(dd->finalize());

	/** set objective bounds */
	primobj_ = dd->getPrimalObjectiveValue();
	dualobj_ = dd->getDualObjectiveValue();
	DSPdebugMessage("primobj %+e, dualobj %+e\n", primobj_, dualobj_);
	message_->print(1, "Best lower bound %e, time elapsed: %.2f sec.\n", dualobj_, dd->getWallTime());

	/** TODO copy primal solution */

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

DSP_RTN_CODE BdDriverMpi::collectSolution()
{
#define FREE_MEMORY            \
	FREE_ARRAY_PTR(nsubprobs)  \
	FREE_ARRAY_PTR(subindices) \
	FREE_ARRAY_PTR(displs)     \
	FREE_ARRAY_PTR(sizesols)   \
	FREE_ARRAY_PTR(sizesolsp)  \
	FREE_ARRAY_PTR(subsols)

	int * nsubprobs  = NULL;
	int * subindices = NULL;
	int * displs     = NULL;
	int * sizesols   = NULL; /**< size of solution vectors for each subproblem */
	int * sizesolsp  = NULL; /**< size of solution vectors for each process */
	double * subsols = NULL;

	BGN_TRY_CATCH

	if (comm_rank_ == 0)
	{
		/** retrieve master pointer */
		BdMaster * master = mw_->getMasterPtr();
		/** collect solutions from the master */
		if (master)
		{
			status_ = master->getStatus();
			primobj_ = master->getPrimalObjective();
			dualobj_ = master->getDualObjective();
			CoinCopyN(master->getPrimalSolution(), master->getSiPtr()->getNumCols(), primsol_);
			numNodes_ = master->getSiPtr()->getNumNodes();
			numIterations_ = master->getSiPtr()->getIterationCount();
		}
		/** nullify master pointer */
		master = NULL;
		DSPdebugMessage("status %d, primobj %+e, dualobj %+e, numNodes %d, numIterations %d\n",
				status_, primobj_, dualobj_, numNodes_, numIterations_);

		/** collect solutions from the workers */

		/** allocate memory */
		nsubprobs  = new int [comm_size_];
		subindices = new int [model_->getNumSubproblems()];
		displs     = new int [comm_size_];
		sizesols   = new int [model_->getNumSubproblems()];
		sizesolsp  = new int [comm_size_];

		/** msg#1. receive number of subproblems */
		int num = 0;
		MPI_Gather(&num, 1, MPI_INT, nsubprobs, 1, MPI_INT, 0, comm_);
		DSPdebugMessage("nsubprobs:\n");
		DSPdebug(DspMessage::printArray(comm_size_, nsubprobs));

		/** set displacement */
		displs[0] = 0;
		for (int i = 1; i < comm_size_; ++i)
			displs[i] = displs[i-1] + nsubprobs[i-1];

		/** msg#2. receive subproblem indices */
		MPI_Gatherv(NULL, 0, MPI_INT, subindices, nsubprobs, displs, MPI_INT, 0, comm_);

		/** msg#3. receive size of solution vectors for each subproblem */
		MPI_Gatherv(NULL, 0, MPI_INT, sizesols, nsubprobs, displs, MPI_INT, 0, comm_);

		/** msg#4. receive size of solution vectors for each process */
		MPI_Gather(&num, 1, MPI_INT, sizesolsp, 1, MPI_INT, 0, comm_);
		DSPdebug(DspMessage::printArray(comm_size_, nsubprobs));

		/** set displacement */
		displs[0] = 0;
		for (int i = 1; i < comm_size_; ++i)
			displs[i] = displs[i-1] + sizesolsp[i-1];
		subsols = new double [displs[comm_size_-1]+sizesolsp[comm_size_-1]];

		/** msg#5. receive subproblem solutions */
		MPI_Gatherv(NULL, 0, MPI_DOUBLE, subsols, sizesolsp, displs, MPI_DOUBLE, 0, comm_);

		/** copy master solution */
		CoinCopyN(mw_->getMasterPtr()->getPrimalSolution(), model_->getNumSubproblemCouplingCols(0), primsol_);

		vector<double*> vsubsols(model_->getNumSubproblems());
		vector<int> vsubsollens(model_->getNumSubproblems());
		for (int i = 0, j = 0; i < model_->getNumSubproblems(); ++i)
		{
			vsubsols[subindices[i]] = subsols + j;
			vsubsollens[subindices[i]] = sizesols[i];
			j += sizesols[i];
		}

		/** copy subproblem solution */
		int pos = model_->getNumSubproblemCouplingCols(0);
		for (int i = 0; i < model_->getNumSubproblems(); ++i)
		{
			CoinCopyN(vsubsols[i], vsubsollens[i], primsol_ + pos);
			pos += vsubsollens[i];
		}
	}
	else
	{
		/** retrieve Benders subproblem pointer */
		BdSub * bdsub = mw_->getWorkerPtr()->getBdSubPtr();
		DSPdebugMessage("bdsub %p\n", bdsub);

		/** msg#1. send number of subproblems */
		int num = bdsub->getNumSubprobs();
		DSPdebugMessage("bdsub->getNumSubprobs() %d\n", num);
		MPI_Gather(&num, 1, MPI_INT, NULL, 0, MPI_INT, 0, comm_);

		/** mgs#2. send subproblem indices */
		MPI_Gatherv(bdsub->getSubprobIndices(), bdsub->getNumSubprobs(),
				MPI_INT, NULL, NULL, NULL, MPI_INT, 0, comm_);

		/** msg#3. send size of solution vectors for each subproblem*/
		int scount = 0;
		sizesols = new int [num];
		for (int i = 0; i < bdsub->getNumSubprobs(); ++i)
		{
			sizesols[i] = bdsub->getNumCols(i);
			scount += bdsub->getNumCols(i);
		}
		MPI_Gatherv(sizesols, bdsub->getNumSubprobs(), MPI_INT, NULL, NULL, NULL, MPI_INT, 0, comm_);

		/** msg#4. send size of solution vectors for each process */
		MPI_Gather(&scount, 1, MPI_INT, NULL, 0, MPI_INT, 0, comm_);

		/** msg#5. send subproblem solutions */
		subsols = new double [scount];
		for (int i = 0, pos = 0; i < bdsub->getNumSubprobs(); ++i)
		{
			CoinCopyN(bdsub->getSolution(i), bdsub->getNumCols(i), subsols + pos);
			pos += bdsub->getNumCols(i);
		}
		MPI_Gatherv(subsols, scount, MPI_DOUBLE, NULL, NULL, NULL, MPI_DOUBLE, 0, comm_);
	}

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}
