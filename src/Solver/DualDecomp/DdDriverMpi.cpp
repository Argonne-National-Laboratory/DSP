/*
 * DdDriverMpi.cpp
 *
 *  Created on: Jul 8, 2016
 *      Author: kibaekkim
 */

#define DSP_DEBUG
#include "DdDriverMpi.h"
#include "Solver/DualDecomp/DdMWSync.h"
#include "Solver/DualDecomp/DdMWAsync.h"

DdDriverMpi::DdDriverMpi(
		DspParams* par,
		DecModel* model,
		MPI_Comm comm):
DdDriver(par, model),
comm_(comm)
{
	BGN_TRY_CATCH

	DSPdebugMessage("comm %d\n", comm_);
	MPI_Comm_rank(comm_, &comm_rank_); /**< get process ID */
	MPI_Comm_size(comm_, &comm_size_); /**< get number of processes */
	DSPdebugMessage("comm_rank_ %d, comm_size_ %d\n", comm_rank_, comm_size_);

	END_TRY_CATCH(;)
}

DSP_RTN_CODE DdDriverMpi::init()
{
	BGN_TRY_CATCH

	/** create Master-Worker framework */
	if (comm_size_ > 1)
	{
		if (par_->getBoolParam("DD/ASYNC"))
			mw_ = new DdMWAsync(comm_, model_, par_, message_);
		else
			mw_ = new DdMWSync(comm_, model_, par_, message_);
		DSP_RTN_CHECK_THROW(mw_->init());
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

		/** primal solution */
		nprimsol = model_->getNumCouplingCols();
		primsol_ = new double [model_->getNumCouplingCols()];
		CoinCopyN(mw_->master_->getBestPrimalSolution(), nprimsol, primsol_);

		/** dual solution */
		ndualsol = model_->getNumCouplingRows();
		dualsol_ = new double [model_->getNumCouplingRows()];
		CoinCopyN(mw_->master_->getBestDualSolution(), ndualsol, dualsol_);

		numNodes_ = mw_->master_->getSiPtr()->getNumNodes();
		numIterations_ = mw_->master_->getSiPtr()->getIterationCount();

		/** communicate */
		MPI_Bcast(&status_, 1, MPI_INT, 0, comm_);
		MPI_Bcast(&primobj_, 1, MPI_DOUBLE, 0, comm_);
		MPI_Bcast(&dualobj_, 1, MPI_DOUBLE, 0, comm_);

		MPI_Bcast(&nprimsol, 1, MPI_INT, 0, comm_);
		MPI_Bcast(primsol_, nprimsol, MPI_DOUBLE, 0, comm_);

		MPI_Bcast(&ndualsol, 1, MPI_INT, 0, comm_);
		MPI_Bcast(dualsol_, ndualsol, MPI_DOUBLE, 0, comm_);

		MPI_Bcast(&numNodes_, 1, MPI_INT, 0, comm_);
		MPI_Bcast(&numIterations_, 1, MPI_INT, 0, comm_);
	}
	else
	{
		MPI_Bcast(&status_, 1, MPI_INT, 0, comm_);
		MPI_Bcast(&primobj_, 1, MPI_DOUBLE, 0, comm_);
		MPI_Bcast(&dualobj_, 1, MPI_DOUBLE, 0, comm_);

		MPI_Bcast(&nprimsol, 1, MPI_INT, 0, comm_);
		primsol_ = new double [nprimsol];
		MPI_Bcast(primsol_, nprimsol, MPI_DOUBLE, 0, comm_);

		MPI_Bcast(&ndualsol, 1, MPI_INT, 0, comm_);
		dualsol_ = new double [ndualsol];
		MPI_Bcast(dualsol_, ndualsol, MPI_DOUBLE, 0, comm_);

		MPI_Bcast(&numNodes_, 1, MPI_INT, 0, comm_);
		MPI_Bcast(&numIterations_, 1, MPI_INT, 0, comm_);
	}

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
