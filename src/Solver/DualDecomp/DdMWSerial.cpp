/*
 * DdMWSerial.cpp
 *
 *  Created on: Apr 8, 2016
 *      Author: kibaekkim
 */

//#define DSP_DEBUG

#include "Solver/DualDecomp/DdMWSerial.h"
#include "Solver/DualDecomp/DdMasterTr.h"
#include "Solver/DualDecomp/DdMasterDsb.h"
#include "Solver/DualDecomp/DdMasterSubgrad.h"

DdMWSerial::DdMWSerial(
		DecModel *   model,  /**< model pointer */
		DspParams *  par,    /**< parameters */
		DspMessage * message /**< message pointer */):
DdMW(model, par, message)
{
	// TODO Auto-generated constructor stub
}

DdMWSerial::~DdMWSerial()
{
	// TODO Auto-generated destructor stub
}

DSP_RTN_CODE DdMWSerial::init()
{
	BGN_TRY_CATCH

	DdMW::init();

	/** set parameters */
	par_->setIntPtrParamSize("ARR_PROC_IDX", model_->getNumSubproblems());
	for (int s = 0; s < model_->getNumSubproblems(); ++s)
		par_->setIntPtrParam("ARR_PROC_IDX", s, s);

	/** create master */
	switch (par_->getIntParam("DD/MASTER_ALGO"))
	{
	case Simplex:
	case IPM:
	case IPM_Feasible:
		master_ = new DdMasterTr(par_, model_, message_);
		break;
	case DSBM:
		master_ = new DdMasterDsb(par_, model_, message_);
		break;
	case Subgradient:
		master_ = new DdMasterSubgrad(par_, model_, message_);
		break;
	}
	DSPdebugMessage("Created master\n");
	/** initialize master */
	master_->init();
	DSPdebugMessage("Initialized master\n");

	/** create LB worker */
	worker_.push_back(new DdWorkerLB(par_, model_, message_));
	/** create CG worker */
	if (parFeasCuts_ >= 0 || parOptCuts_ >= 0)
		worker_.push_back(new DdWorkerCGBd(par_, model_, message_));
	/** create UB worker */
	if (parEvalUb_ >= 0)
		worker_.push_back(new DdWorkerUB(par_, model_, message_));
	/** initialize workers */
	for (unsigned i = 0; i < worker_.size(); ++i)
		worker_[i]->init();

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMWSerial::finalize()
{
	BGN_TRY_CATCH
#if 0
	char filename[64];
	const char * output_prefix = par_->getStrParam("OUTPUT/PREFIX").c_str();

	/** free master */
	sprintf(filename, "%s-DD.out", output_prefix);
	writeIterInfo(filename);
	sprintf(filename, "%s-Master.out", output_prefix);
	master_->write(filename);
#endif
	master_->finalize();
	FREE_PTR(master_);

	/** free workers */
	for (unsigned i = 0; i < worker_.size(); ++i)
	{
#if 0
		switch (worker_[i]->getType())
		{
		case DdWorker::LB:
			sprintf(filename, "%s-LB.out", output_prefix);
			break;
		case DdWorker::UB:
			sprintf(filename, "%s-UB.out", output_prefix);
			break;
		case DdWorker::CGBd:
			sprintf(filename, "%s-CG.out", output_prefix);
			break;
		}
		worker_[i]->write(filename);
#endif
		worker_[i]->finalize();
		FREE_PTR(worker_[i]);
	}
	worker_.clear();

	/** finalize master-worker */
	DdMW::finalize();

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMWSerial::run()
{
#define FREE_MEMORY \
	FREE_ARRAY_PTR(lambdas) \
	thetas = NULL;

	const double * thetas  = NULL; /**< of master problem */
	double **      lambdas = NULL; /**< of master problem */

	/** pointers to workers */
	DdWorkerLB * workerlb = NULL;
	DdWorkerCGBd * workercg = NULL;
	DdWorkerUB * workerub = NULL;

	Solutions coupling_solutions; /**< coupling solutions */

	/** Benders cuts */
	int ncuts = 0;
	OsiCuts cuts, emptycuts;
	int cg_status = DSP_STAT_MW_CONTINUE;

	BGN_TRY_CATCH

	/** retrieve DdWorker pointers */
	for (unsigned i = 0; i < worker_.size(); ++i)
	{
		switch(worker_[i]->getType())
		{
		case DdWorker::LB:
			workerlb = dynamic_cast<DdWorkerLB*>(worker_[i]);
			break;
		case DdWorker::CGBd:
			workercg = dynamic_cast<DdWorkerCGBd*>(worker_[i]);
			break;
		case DdWorker::UB:
			workerub = dynamic_cast<DdWorkerUB*>(worker_[i]);
			break;
		default:
			message_->print(0, "Unknown worker type (%d).\n", worker_[i]->getType());
			break;
		}
	}
	assert(workerlb!=NULL);

	/** allocate memory for lambdas */
	lambdas = new double * [model_->getNumSubproblems()];

	printHeaderInfo();

	/** reset iteration info */
	itercnt_   = 0;
	iterstime_ = CoinGetTimeOfDay();

	/**
	 * This is the main loop to iteratively solve master problem.
	 */
	while (1)
	{
		/** reset cut generation status */
		cg_status = DSP_STAT_MW_CONTINUE;

		/** reset iteration code */
		itercode_ = ' ';

		/** Solve subproblems */
		DSP_RTN_CHECK_THROW(workerlb->solve());

		/** update master */
		double subprimobj = 0.0;
		double subdualobj = 0.0;
		for (int s = 0; s < model_->getNumSubproblems(); ++s)
		{
			int sindex = workerlb->subprobs_[s]->sind_;
			master_->subprimobj_[sindex] = workerlb->subprobs_[s]->getPrimalBound();
			master_->subdualobj_[sindex] = workerlb->subprobs_[s]->getDualBound();
			CoinCopyN(workerlb->subprobs_[s]->si_->getSolution(),
					workerlb->subprobs_[s]->ncols_coupling_,
					master_->subsolution_[sindex]);
			subprimobj += workerlb->subprobs_[s]->getPrimalBound();
			subdualobj += workerlb->subprobs_[s]->getDualBound();
			/*DSPdebugMessage("Scenario %d, primobj %+e, dualobj %+e, primobj_sum %+e, dualobj_sum %+e\n",
					workerlb->subprobs_[s]->sind_,
					workerlb->subprobs_[s]->getPrimalBound(),
					workerlb->subprobs_[s]->getDualBound(),
					subprimobj, subdualobj);*/
		}

		if (parEvalUb_ >= 0 || parFeasCuts_ >= 0 || parOptCuts_ >= 0)
		{
			/** store coupling solutions */
			DSP_RTN_CHECK_THROW(storeCouplingSolutions(coupling_solutions));
			/** generate cuts */
			if (parFeasCuts_ >= 0 || parOptCuts_ >= 0)
			{
				cg_status = generateBendersCuts(workercg, coupling_solutions, cuts);
				/** resolve subproblems? */
				if (cg_status == DSP_STAT_MW_RESOLVE)
				{
					printIterInfo();
					DSPdebugMessage("Resolve subproblems.\n");
					continue;
				}
				/** move cuts to a global pool */
				for (int i = 0; i < cuts.sizeCuts(); ++i)
				{
					OsiRowCut * rc = cuts.rowCutPtr(i);
					cutsToAdd_->insert(rc);
				}
				cuts.dumpCuts();
			}
			/** evaluate coupling solutions */
			if (parEvalUb_ >= 0)
			{
				double oldub = master_->bestprimobj_;
				for (unsigned i = 0; i < coupling_solutions.size(); ++i)
				{
					double newub = workerub->evaluate(coupling_solutions[i]);
					DSPdebugMessage("Current upper bound %+e\n", newub);
					if (newub < master_->bestprimobj_)
						master_->bestprimobj_ = newub;
				}
				if (oldub > master_->bestprimobj_)
					itercode_ = 'P';
			}
			/** clear stored solutions */
			coupling_solutions.clear();
		}

		/** update problem */
		double olddual = master_->bestdualobj_;
		master_->updateProblem();
		if (olddual < master_->bestdualobj_)
			itercode_ = itercode_ == 'P' ? 'B' : 'D';

		/** STOP with small gap */
		if (itercnt_ >= master_->getParPtr()->getIntParam("DD/ITER_LIM"))
		{
			message_->print(1, "The iteration limit is reached.\n");
			break;
		}

		/** solve problem */
		double tic = CoinGetTimeOfDay();
		DSP_RTN_CHECK_THROW(master_->solve());
		DSPdebugMessage("Solved the master (%.2f sec).\n", CoinGetTimeOfDay() - tic);

		printIterInfo();

		/** termination test */
		if (master_->terminationTest() == DSP_STAT_MW_STOP) break;
		/** check duality gap */
		if (master_->getRelDualityGap() < master_->getParPtr()->getDblParam("DD/STOP_TOL"))
		{
			message_->print(1, "The duality gap limit %+e is reached.\n", master_->getRelDualityGap());
			break;
		}

		/** increment iteration count */
		itercnt_++;

		/** retrieve master solution by part */
		double * master_primsol = const_cast<double*>(master_->getPrimalSolution());
		thetas  = master_primsol;
		for (int i = 0, j = model_->getNumSubproblems(); i < model_->getNumSubproblems(); ++i)
		{
			/** shallow copy */
			lambdas[i] = master_primsol + j;
			j += model_->getNumSubproblemCouplingRows(i);
		}
		master_primsol = NULL;

		/** update subproblems */
		for (int s = 0; s < model_->getNumSubproblems(); ++s)
		{
			int sindex = workerlb->subprobs_[s]->sind_;
			workerlb->subprobs_[s]->theta_ = thetas[sindex];
			workerlb->subprobs_[s]->updateProblem(lambdas[sindex], master_->bestprimobj_);
			/** apply Benders cuts */
			if (parFeasCuts_ >= 0 || parOptCuts_ >= 0)
			{
				workerlb->subprobs_[s]->pushCuts(cutsToAdd_);
				DSPdebugMessage("Pushed %d Benders cuts.\n", cutsToAdd_->sizeCuts());
			}
		}
	}

	/** set best dual objective */
	//master_->bestdualobj_ = master_->primobj_;

	DSPdebugMessage2("primsol_:\n");
	DSPdebug2(message_->printArray(master_->getSiPtr()->getNumCols(), master_->getPrimalSolution()));

	/** release shallow-copy of pointers */
	for (int i = 0; i < model_->getNumSubproblems(); ++i)
		lambdas[i] = NULL;

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

DSP_RTN_CODE DdMWSerial::generateBendersCuts(
		DdWorkerCGBd * workercg, /**< CG worker pointer */
		Solutions solutions, /**< solutions at which cuts are generated */
		OsiCuts& cuts)
{
#define FREE_MEMORY \
	FREE_ARRAY_PTR(aggcut)

	int ret = DSP_STAT_MW_CONTINUE;

	if (solutions.size() == 0) return ret;
	if (model_->isStochastic() == false)
	{
		message_->print(0, "This problem is not a stochastic program. Benders cut option is valid only for stochastic programming.\n");
		parFeasCuts_ = -1;
		parOptCuts_ = -1;
		return ret;
	}

	vector<int> cuttype;
	TssModel * tssmodel = NULL;
	double * aggcut = NULL;
	double aggrhs;

	BGN_TRY_CATCH

	/** timing */
	double cputime = CoinCpuTime();
	double walltime = CoinGetTimeOfDay();

	/** downcast to TssModel */
	tssmodel = dynamic_cast<TssModel*>(model_);

	/** allocate memory */
	aggcut = new double [tssmodel->getNumCols(0) + 1];

	/** generate cuts */
	for (unsigned i = 0; i < solutions.size(); ++i)
	{
		/** create local cut pools */
		OsiCuts localcuts;

		/** generate cuts at a solution */
		DSP_RTN_CHECK_THROW(workercg->generateCuts(solutions[i], &localcuts, cuttype));
		DSPdebugMessage("Benders cut generator generated %d cuts for solution %d.\n",
				localcuts.sizeCuts(), i);

		/** check if there exists a feasibility cut */
		bool hasFeasibilityCuts = false;
		bool hasOptimalityCuts = true;
		for (unsigned j = 0; j < cuttype.size(); ++j)
		{
			if (cuttype[j] == DdWorkerCGBd::Feas)
			{
				hasFeasibilityCuts = true;
				hasOptimalityCuts = false;
				break;
			}
			if (cuttype[j] != DdWorkerCGBd::Opt)
				hasOptimalityCuts = false;
		}

		if ((hasFeasibilityCuts == true && parFeasCuts_ >= 0) ||
			(hasOptimalityCuts == true && parOptCuts_ >= 0))
		{
			if (hasOptimalityCuts)
			{
				/** construct optimality cut */
				CoinZeroN(aggcut, tssmodel->getNumCols(0) + 1);
				for (int j = 0; j < tssmodel->getNumCols(0); ++j)
					aggcut[j] = -(tssmodel->getObjCore(0)[j]);
				aggrhs = 0.0;

				for (int j = 0; j < localcuts.sizeCuts(); ++j)
				{
					OsiRowCut * rc = localcuts.rowCutPtr(j);
					CoinPackedVector cutrow = rc->row();
					for (int k = 0; k < cutrow.getNumElements(); ++k)
						aggcut[cutrow.getIndices()[k]] += cutrow.getElements()[k];
					aggrhs += rc->lb();
					DSPdebugMessage2("localcuts[%d]:\n", j);
					DSPdebug2(rc->print());
				}
				DSPdebug2(message_->printArray(tssmodel->getNumCols(0) + 1, aggcut));

				CoinPackedVector orow;
				for (int j = 0; j < tssmodel->getNumCols(0); ++j)
					if (fabs(aggcut[j]) > 1.0e-8)
						orow.insert(j, aggcut[j]);
				orow.insert(tssmodel->getNumCols(0) + tssmodel->getNumCols(1), 1.0);

				OsiRowCut ocut;
				ocut.setRow(orow);
				ocut.setUb(COIN_DBL_MAX);
				ocut.setLb(aggrhs);

				DSPdebug2(ocut.print());
				cuts.insertIfNotDuplicate(ocut);
			}

			if (hasFeasibilityCuts)
			{
				/** copy cut pointers */
				for (int j = 0; j < localcuts.sizeCuts(); ++j)
				{
					OsiRowCut * rc = localcuts.rowCutPtr(j);
					DSPdebug(rc->print());
					cuts.insertIfNotDuplicate(*rc);
					ret = DSP_STAT_MW_RESOLVE;
				}
			}
		}

		/** free cuts */
		for (int j = 0; j < localcuts.sizeCuts(); ++j)
		{
			OsiRowCut * rc = localcuts.rowCutPtr(j);
			FREE_PTR(rc);
		}
		localcuts.dumpCuts();
	}

	/** log timing results */
	workercg->s_statuses_.push_back(DSP_RTN_OK);
	workercg->s_primobjs_.push_back(0.0);
	workercg->s_dualobjs_.push_back(0.0);
	workercg->s_cputimes_.push_back(CoinCpuTime() - cputime);
	workercg->s_walltimes_.push_back(CoinGetTimeOfDay() - walltime);

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return ret;
#undef FREE_MEMORY
}
