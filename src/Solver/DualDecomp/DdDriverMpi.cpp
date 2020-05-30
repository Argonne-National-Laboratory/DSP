/*
 * DdDriverMpi.cpp
 *
 *  Created on: Jul 8, 2016
 *      Author: kibaekkim
 */

// #define DSP_DEBUG
#include "DdDriverMpi.h"
#include "Solver/DualDecomp/DdMWSync.h"
#include "Solver/DualDecomp/DdMWAsync.h"
#include "Solver/DualDecomp/DdMWAsyncDyn.h"

DdDriverMpi::DdDriverMpi(
			DecModel*   model,   /**< model pointer */
			DspParams*  par,     /**< parameters */
			DspMessage* message, /**< message pointer */
			MPI_Comm    comm     /**< MPI communicator */):
DdDriver(model, par, message),
comm_(comm)
{
	BGN_TRY_CATCH

	MPI_Comm_rank(comm_, &comm_rank_); /**< get process ID */
	MPI_Comm_size(comm_, &comm_size_); /**< get number of processes */
	DSPdebugMessage("comm_rank_ %d, comm_size_ %d\n", comm_rank_, comm_size_);

	END_TRY_CATCH(;)
}

DdDriverMpi::~DdDriverMpi() {}

DdDriverMpi::DdDriverMpi(const DdDriverMpi& rhs) :
DdDriver(rhs),
comm_(rhs.comm_),
comm_rank_(rhs.comm_rank_),
comm_size_(rhs.comm_size_) {}

DSP_RTN_CODE DdDriverMpi::init()
{
	BGN_TRY_CATCH

	/** create Master-Worker framework */
	if (comm_size_ > 1)
	{
		if (par_->getBoolParam("DD/ASYNC")) {
			if (par_->getBoolParam("DD/ASYNC/DYNAMIC"))
				mw_ = new DdMWAsyncDyn(comm_, model_, par_, message_);
			else
				mw_ = new DdMWAsync(comm_, model_, par_, message_);
		} else
			mw_ = new DdMWSync(comm_, model_, par_, message_);
		DSP_RTN_CHECK_THROW(mw_->init());
	} else {
		printf("Number of processes should be greater than 1.\n");
		return DSP_RTN_ERR;
	}

	DdDriver::init();

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdDriverMpi::run()
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

	/** collect solutions */
	int nprimsol, ndualsol;
	if (comm_rank_ == 0)
	{
		status_ = mw_->master_->getStatus();
		primobj_ = mw_->master_->getBestPrimalObjective();
		dualobj_ = mw_->master_->getBestDualObjective();
		bestprimobj_ = primobj_;
		bestdualobj_ = dualobj_;

		/** primal solution */
		nprimsol = model_->getFullModelNumCols();
		primsol_.resize(nprimsol);
		CoinCopyN(mw_->master_->getBestPrimalSolution(), nprimsol, &primsol_[0]);

		/** dual solution */
		ndualsol = model_->getNumCouplingRows();
		dualsol_.resize(ndualsol);
		CoinCopyN(mw_->master_->getBestDualSolution(), ndualsol, &dualsol_[0]);

		numNodes_ = mw_->master_->getDspOsiPtr()->getNumNodes();
		numIterations_ = mw_->master_->getSiPtr()->getIterationCount();

		/** communicate */
		MPI_Bcast(&status_, 1, MPI_INT, 0, comm_);
		MPI_Bcast(&primobj_, 1, MPI_DOUBLE, 0, comm_);
		MPI_Bcast(&dualobj_, 1, MPI_DOUBLE, 0, comm_);

		MPI_Bcast(&nprimsol, 1, MPI_INT, 0, comm_);
		MPI_Bcast(&primsol_[0], nprimsol, MPI_DOUBLE, 0, comm_);

		MPI_Bcast(&ndualsol, 1, MPI_INT, 0, comm_);
		MPI_Bcast(&dualsol_[0], ndualsol, MPI_DOUBLE, 0, comm_);

		MPI_Bcast(&numNodes_, 1, MPI_INT, 0, comm_);
		MPI_Bcast(&numIterations_, 1, MPI_INT, 0, comm_);
	}
	else
	{
		MPI_Bcast(&status_, 1, MPI_INT, 0, comm_);
		MPI_Bcast(&primobj_, 1, MPI_DOUBLE, 0, comm_);
		MPI_Bcast(&dualobj_, 1, MPI_DOUBLE, 0, comm_);

		MPI_Bcast(&nprimsol, 1, MPI_INT, 0, comm_);
		primsol_.resize(nprimsol);
		MPI_Bcast(&primsol_[0], nprimsol, MPI_DOUBLE, 0, comm_);

		MPI_Bcast(&ndualsol, 1, MPI_INT, 0, comm_);
		dualsol_.resize(ndualsol);
		MPI_Bcast(&dualsol_[0], ndualsol, MPI_DOUBLE, 0, comm_);

		MPI_Bcast(&numNodes_, 1, MPI_INT, 0, comm_);
		MPI_Bcast(&numIterations_, 1, MPI_INT, 0, comm_);
	}
	bestprimsol_ = primsol_;
	bestdualsol_ = dualsol_;

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdDriverMpi::finalize()
{
	BGN_TRY_CATCH

	/** finalize master-worker framework */
	mw_->finalize();

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}
