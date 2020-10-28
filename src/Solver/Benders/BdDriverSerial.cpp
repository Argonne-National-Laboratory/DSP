/*
 * BdDriverSerial.cpp
 *
 *  Created on: Jul 10, 2016
 *      Author: kibaekkim
 */

// #define DSP_DEBUG

#include "Solver/Benders/BdDriverSerial.h"
#include "Solver/Benders/BdMWSerial.h"
#include "Solver/DualDecomp/DdDriverSerial.h"

BdDriverSerial::BdDriverSerial(
			DecModel *   model,  /**< model pointer */
			DspParams *  par,    /**< parameters */
			DspMessage * message /**< message pointer */):
BdDriver(model,par,message) {}

BdDriverSerial::BdDriverSerial(const BdDriverSerial& rhs) :
BdDriver(rhs) {}

BdDriverSerial::~BdDriverSerial() {}

DSP_RTN_CODE BdDriverSerial::init()
{
	BGN_TRY_CATCH

	/** primal soltuion */
	primsol_.resize(model_->getFullModelNumCols());

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
		DSP_RTN_CHECK_THROW(findLowerBound());

	/** set objective bounds */
	DSPdebugMessage("setObjectiveBounds\n");
	mw_->getMasterPtr()->setObjectiveBounds(primobj_, dualobj_);

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
	dd = new DdDriverSerial(model_, par_, message_);
	DSP_RTN_CHECK_THROW(dd->init());
	DSP_RTN_CHECK_THROW(dd->run());
	DSP_RTN_CHECK_THROW(dd->finalize());

	/** set objective bounds */
	primobj_ = dd->getPrimalObjective();
	dualobj_ = dd->getDualObjective();
	DSPdebugMessage("primobj %+e, dualobj %+e\n", primobj_, dualobj_);
	message_->print(1, "Best lower bound %e, time elapsed: %.2f sec.\n", dualobj_, dd->getWallTime());

	/** TODO copy primal solution */

	/** rollback parameters */
	par_->setIntParam("DD/ITER_LIM", iterlim);
	par_->setIntParam("DD/FEAS_CUTS", fcut);
	par_->setIntParam("DD/OPT_CUTS", ocut);
	par_->setIntParam("DD/EVAL_UB", evalub);
	model_->setDro(isdro);

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
		bestprimobj_ = primobj_;
		dualobj_ = master->getDualObjective();
		bestdualobj_ = dualobj_;
		DSPdebugMessage("status %d, primobj %+e, dualobj %+e\n", status_, primobj_, dualobj_);
		CoinCopyN(master->getPrimalSolution(), model_->getNumSubproblemCouplingCols(0), &primsol_[0]);
		bestprimsol_ = primsol_;
		numNodes_ = master->getNumNodes();
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
		CoinCopyN(bdsub->getSolution(i), bdsub->getNumCols(i), &primsol_[pos]);
		pos += bdsub->getNumCols(i);
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}
