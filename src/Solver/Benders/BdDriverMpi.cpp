/*
 * BdDriverMpi.cpp
 *
 *  Created on: Jul 10, 2016
 *      Author: kibaekkim
 */

//#define DSP_DEBUG
#include "Solver/Benders/BdDriverMpi.h"
#include "Solver/DualDecomp/DdDriverMpi.h"

BdDriverMpi::BdDriverMpi(
		DecModel *   model,   /**< model pointer */
		DspParams *  par,     /**< parameters */
		DspMessage * message, /**< message pointer */
		MPI_Comm comm)        /**< communicator */:
BdDriver(model, par, message),
comm_(comm) {
	MPI_Comm_rank(comm, &comm_rank_);
	MPI_Comm_size(comm, &comm_size_);
}

BdDriverMpi::BdDriverMpi(const BdDriverMpi& rhs) :
BdDriver(rhs),
comm_(rhs.comm_),
comm_rank_(rhs.comm_rank_),
comm_size_(rhs.comm_size_) {}

BdDriverMpi::~BdDriverMpi() {}

DSP_RTN_CODE BdDriverMpi::init()
{
	BGN_TRY_CATCH

	/** allocate memory for primal solution */
	primsol_.resize(model_->getFullModelNumCols());

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
		/** set objective bounds */
		DSPdebugMessage("setObjectiveBounds\n");
		mw_->getMasterPtr()->setObjectiveBounds(primobj_, dualobj_);

		/** set initial solutions */
		if (initsols_.size() > 0)
			DSP_RTN_CHECK_THROW(mw_->getMasterPtr()->setSolutions(initsols_));
        /** set time limit */
        mw_->getMasterPtr()->setTimeLimit(par_->getDblParam("BD/WALL_LIM") - CoinGetTimeOfDay() + walltime_);
	}

	/** run */
	DSP_RTN_CHECK_THROW(mw_->run());

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

DSP_RTN_CODE BdDriverMpi::findLowerBound() {
#define FREE_MEMORY \
    FREE_PTR(dd);

    DdDriver *dd = NULL;

    BGN_TRY_CATCH

	/** set parameters */
	int iterlim = par_->getIntParam("DD/ITER_LIM");
	int fcut = par_->getIntParam("DD/FEAS_CUTS");
	int ocut = par_->getIntParam("DD/OPT_CUTS");
	int evalub = par_->getIntParam("DD/EVAL_UB");
	int iter_lim = model_->isDro() ? 0 : par_->getIntParam("BD/DD/ITER_LIM");
	bool isdro = model_->isDro();
	par_->setIntParam("DD/ITER_LIM", iter_lim);
	par_->setIntParam("DD/FEAS_CUTS", -1);
	par_->setIntParam("DD/OPT_CUTS", -1);
	par_->setIntParam("DD/EVAL_UB", -1);
	model_->setDro(false);

	message_->print(1, "Finding a good lower bound using Dual Decomposition...\n");

	/** use dual decomposition */
	dd = new DdDriverMpi(model_, par_, message_, comm_);
	DSP_RTN_CHECK_THROW(dd->init());
	DSP_RTN_CHECK_THROW(dd->run());
	DSP_RTN_CHECK_THROW(dd->finalize());

	/** set objective bounds */
	primobj_ = dd->getPrimalObjective();
	dualobj_ = dd->getDualObjective();
	message_->print(1, "Best lower bound %e, time elapsed: %.2f sec.\n", dualobj_, dd->getWallTime());
	DSPdebugMessage("Rank %d: primobj %+e, dualobj %+e\n", comm_rank_, primobj_, dualobj_);

	/** TODO copy primal solution */

	/** rollback parameters */
	par_->setIntParam("DD/ITER_LIM", iterlim);
	par_->setIntParam("DD/FEAS_CUTS", fcut);
	par_->setIntParam("DD/OPT_CUTS", ocut);
	par_->setIntParam("DD/EVAL_UB", evalub);
	model_->setDro(isdro);

	END_TRY_CATCH_RTN(FREE_MEMORY, DSP_RTN_ERR)

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
		/** collect solutions from the master */
        status_ = mw_->getMasterPtr()->getStatus();
        primobj_ = mw_->getMasterPtr()->getPrimalObjective();
        dualobj_ = mw_->getMasterPtr()->getDualObjective();
        numNodes_ = mw_->getMasterPtr()->getNumNodes();
        numIterations_ = mw_->getMasterPtr()->getSiPtr()->getIterationCount();
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
        DSPdebugMessage("sizesolsp:\n");
		DSPdebug(DspMessage::printArray(comm_size_, sizesolsp));

		/** set displacement */
		displs[0] = 0;
		for (int i = 1; i < comm_size_; ++i)
			displs[i] = displs[i-1] + sizesolsp[i-1];
		subsols = new double [displs[comm_size_-1]+sizesolsp[comm_size_-1]];

		/** msg#5. receive subproblem solutions */
		MPI_Gatherv(NULL, 0, MPI_DOUBLE, subsols, sizesolsp, displs, MPI_DOUBLE, 0, comm_);

		/** copy master solution */
		CoinCopyN(mw_->getMasterPtr()->getPrimalSolution(), model_->getNumSubproblemCouplingCols(0), &primsol_[0]);

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
			CoinCopyN(vsubsols[subindices[i]], vsubsollens[subindices[i]], &primsol_[pos]);
			pos += vsubsollens[subindices[i]];
		}
        DSPdebugMessage("primsol_(%d):\n", model_->getFullModelNumCols());
        DSPdebug(DspMessage::printArray(model_->getFullModelNumCols(), &primsol_[0]));

        /** broadcast solutions */
        MPI_Bcast(&status_, 1, MPI_INT, 0, comm_);
        MPI_Bcast(&primobj_, 1, MPI_DOUBLE, 0, comm_);
        MPI_Bcast(&dualobj_, 1, MPI_DOUBLE, 0, comm_);
        MPI_Bcast(&primsol_[0], model_->getFullModelNumCols(), MPI_DOUBLE, 0, comm_);
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
        DSPdebugMessage("subproblem solutions:\n");
        DSPdebug(DspMessage::printArray(scount, subsols));

        /** broadcast solutions */
        MPI_Bcast(&status_, 1, MPI_INT, 0, comm_);
        MPI_Bcast(&primobj_, 1, MPI_DOUBLE, 0, comm_);
        MPI_Bcast(&dualobj_, 1, MPI_DOUBLE, 0, comm_);
        MPI_Bcast(&primsol_[0], model_->getFullModelNumCols(), MPI_DOUBLE, 0, comm_);
	}
	bestprimobj_ = primobj_;
	bestdualobj_ = dualobj_;
	bestprimsol_ = primsol_;

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}
