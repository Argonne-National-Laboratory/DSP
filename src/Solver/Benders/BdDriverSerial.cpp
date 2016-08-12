/*
 * BdDriverSerial.cpp
 *
 *  Created on: Jul 10, 2016
 *      Author: kibaekkim
 */

//#define DSP_DEBUG

#include "Solver/Benders/BdDriverSerial.h"
#include "Solver/Benders/BdMWSerial.h"
#include "Solver/DualDecomp/DdDriverSerial.h"

BdDriverSerial::BdDriverSerial(
		DspParams * par, /**< parameter pointer */
		DecModel * model /**< model pointer */):
BdDriver(par,model)
{
}

BdDriverSerial::~BdDriverSerial()
{
	// TODO Auto-generated destructor stub
}

DSP_RTN_CODE BdDriverSerial::init()
{
	BGN_TRY_CATCH

	/** primal soltuion */
	primsol_  = new double [model_->getFullModelNumCols()];

	/** create Master-Worker framework */
	mw_ = new BdMWSerial(model_, par_, message_);
	DSPdebugMessage("Created BdMWSerial()\n");

	/** initialize master-worker framework */
	mw_->init();
	DSPdebugMessage("Initialized BdMWSerial()\n");

	/** set auxiliary data */
	mw_->getMasterPtr()->setAuxVarData(aux_size_, aux_obj_, aux_clbd_, aux_cubd_);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE BdDriverSerial::run()
{
	BGN_TRY_CATCH

	/** tic */
	cputime_  = CoinCpuTime();
	walltime_ = CoinGetTimeOfDay();

	/** find lower bound first if it is not given by user */
	if (dualobj_ <= -COIN_DBL_MAX)
		findLowerBound();

	/** set objective bounds */
	DSPdebugMessage("setObjectiveBounds\n");
	mw_->getMasterPtr()->setObjectiveBounds(primobj_, dualobj_);

	/** set branching priorities */
	DSPdebugMessage("setBranchingPriority\n");
	if (numPriorities_ > 0)
		mw_->getMasterPtr()->setBranchingPriority(numPriorities_, priorities_);

	/** set initial solutions */
	DSPdebugMessage("setSolutions\n");
	if (initsols_.size() > 0)
		mw_->getMasterPtr()->setSolutions(initsols_);

	/** run */
	mw_->getMasterPtr()->setTimeLimit(
			par_->getDblParam("BD/WALL_LIM") - CoinGetTimeOfDay() + walltime_);
	DSPdebugMessage("Solve Benders\n");
	mw_->run();

	/** toc */
	cputime_  = CoinCpuTime() - cputime_;
	walltime_ = CoinGetTimeOfDay() - walltime_;
	DSPdebugMessage("Time: cpu %.2f wall %.2f\n", cputime_, walltime_);

	/** collect solutions from the master and the worker */
	collectSolution();
	DSPdebugMessage("Collected solutions (status %d).\n", status_);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE BdDriverSerial::finalize()
{
	BGN_TRY_CATCH

	/** finalize master-worker framework */
	mw_->finalize();

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE BdDriverSerial::findLowerBound()
{
#define FREE_MEMORY \
	FREE_PTR(dd);

	DdDriver * dd = NULL;

	BGN_TRY_CATCH

	message_->print(1, "Finding a good lower bound using Dual Decomposition...\n");

	/** use dual decomposition */
	dd = new DdDriverSerial(par_, model_);
	dd->init();
	dd->run();
	dd->finalize();

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

DSP_RTN_CODE BdDriverSerial::collectSolution()
{
	BGN_TRY_CATCH

	/** retrieve master pointer */
	BdMaster * master = mw_->getMasterPtr();
	/** collect solution from the master */
	if (master)
	{
		status_ = master->getStatus();
		primobj_ = master->getPrimalObjective();
		dualobj_ = master->getDualObjective();
		DSPdebugMessage("status %d, primobj %+e, dualobj %+e\n", status_, primobj_, dualobj_);
		CoinCopyN(master->getPrimalSolution(), model_->getNumSubproblemCouplingCols(0), primsol_);
		numNodes_ = master->getSiPtr()->getNumNodes();
		numIterations_ = master->getSiPtr()->getIterationCount();
		DSPdebugMessage("nodes %d, iterations %d\n", numNodes_, numIterations_);
	}
	/** nullify master pointer */
	master = NULL;

	/** retrieve Benders subproblem pointer */
	BdSub * bdsub = mw_->getWorkerPtr()->getBdSubPtr();
	/** copy subproblem solution */
	int pos = model_->getNumSubproblemCouplingCols(0);
	for (int i = 0; i < model_->getNumSubproblems(); ++i)
	{
		//DSPdebugMessage("Subprob %d numcols %d\n", i, bdsub->getNumCols(i));
		CoinCopyN(bdsub->getSolution(i), bdsub->getNumCols(i), primsol_ + pos);
		pos += bdsub->getNumCols(i);
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}
