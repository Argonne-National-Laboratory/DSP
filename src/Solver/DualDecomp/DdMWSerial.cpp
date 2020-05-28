/*
 * DdMWSerial.cpp
 *
 *  Created on: Apr 8, 2016
 *      Author: kibaekkim
 */

// #define DSP_DEBUG

#include "Model/TssModel.h"
#include "Solver/DualDecomp/DdMWSerial.h"
#include "Solver/DualDecomp/DdMasterTr.h"
// #ifdef DSP_HAS_OOQP
// #include "Solver/DualDecomp/DdMasterDsb.h"
// #endif
#include "Solver/DualDecomp/DdMasterSubgrad.h"

DdMWSerial::DdMWSerial(
		DecModel *   model,  /**< model pointer */
		DspParams *  par,    /**< parameters */
		DspMessage * message /**< message pointer */):
DdMW(model, par, message) {}

DdMWSerial::DdMWSerial(const DdMWSerial& rhs) :
DdMW(rhs) {}

DdMWSerial::~DdMWSerial() {}

DSP_RTN_CODE DdMWSerial::init() {
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
		master_ = new DdMasterTr(model_, par_, message_);
		break;
	case DSBM:
		// master_ = new DdMasterDsb(model_, par_, message_);
		char msg[128];
		sprintf(msg, "DD/MASTER_ALGO = %d is not currently supported.", DSBM);
		throw CoinError(msg, "init", "DdMWSerial");
		break;
	case Subgradient:
		master_ = new DdMasterSubgrad(model_, par_, message_);
		break;
	}
	DSPdebugMessage("Created master\n");
	/** initialize master */
	DSP_RTN_CHECK_THROW(master_->init());
	DSPdebugMessage("Initialized master\n");

	/** create LB worker */
	worker_.push_back(new DdWorkerLB(model_, par_, message_));
	/** create CG worker */
#ifdef DSP_HAS_SCIP
	if (par_->getIntParam("DW/SUB/SOLVER") == OsiScip && (parFeasCuts_ >= 0 || parOptCuts_ >= 0))
		worker_.push_back(new DdWorkerCGBd(model_, par_, message_));
	else {
		message_->print(0, "No Benders cut is generated.\n");
		parFeasCuts_ = -1;
		parOptCuts_ = -1;
	}
#endif
	/** create UB worker */
	if (parEvalUb_ >= 0)
		worker_.push_back(new DdWorkerUB(model_, par_, message_));
	/** initialize workers */
	for (unsigned i = 0; i < worker_.size(); ++i)
		worker_[i]->init();

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMWSerial::finalize() {
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

DSP_RTN_CODE DdMWSerial::run() {
#define FREE_MEMORY \
	FREE_ARRAY_PTR(lambdas) \
	FREE_ARRAY_PTR(nsubsolution) \
	thetas = NULL; \
	Ps = NULL;

	const double * thetas  = NULL; /**< of master problem */
	double **      lambdas = NULL; /**< of master problem */
	const double * Ps      = NULL; /**< of DRO master problem */

	int * nsubsolution = NULL; /**< size of subproblem solution */

	/** pointers to workers */
	DdWorkerLB * workerlb = NULL;
#ifdef DSP_HAS_SCIP
	DdWorkerCGBd * workercg = NULL;
#endif
	DdWorkerUB * workerub = NULL;
	TssModel* tss = NULL;

	Solutions coupling_solutions; /**< coupling solutions */

	/** Benders cuts */
	OsiCuts cuts, emptycuts;
	int cg_status = DSP_STAT_MW_CONTINUE;

	BGN_TRY_CATCH

	if (model_->isStochastic()) {
		try {
			tss = dynamic_cast<TssModel*>(model_);
		} catch (const std::bad_cast &e) {
			printf("Error: Model claims to be stochastic when it is not\n");
            return DSP_RTN_ERR;
		}
	}

	/** retrieve DdWorker pointers */
	for (unsigned i = 0; i < worker_.size(); ++i)
	{
		switch(worker_[i]->getType())
		{
		case DdWorker::LB:
			workerlb = dynamic_cast<DdWorkerLB*>(worker_[i]);
			break;
		case DdWorker::CGBd:
#ifdef DSP_HAS_SCIP
			workercg = dynamic_cast<DdWorkerCGBd*>(worker_[i]);
#endif
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
	/** NOTE: lambdas will be used as shallow pointers. */
	lambdas = new double * [model_->getNumSubproblems()];

	nsubsolution = new int [model_->getNumSubproblems()];
	for (int s = 0; s < model_->getNumSubproblems(); ++s) {
		nsubsolution[s] = workerlb->subprobs_[s]->getNumCols();
		DSPdebugMessage("nsubsolution[%d] = %d\n", s, nsubsolution[s]);
	}

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

		/** set time limit */
		workerlb->setTimeLimit(remainingTime());

		/** Solve subproblems */
		DSP_RTN_CHECK_THROW(workerlb->solve());

		/** update master */
		double subprimobj = 0.0;
		double subdualobj = 0.0;
		for (int s = 0; s < model_->getNumSubproblems(); ++s)
		{
			int sindex = workerlb->subprobs_[s]->sind_;
			master_->subprimobj_[sindex] = workerlb->subprobs_[s]->getPrimalObjective();
			master_->subdualobj_[sindex] = workerlb->subprobs_[s]->getDualObjective();
			CoinCopyN(workerlb->subprobs_[s]->getPrimalSolution(),
					nsubsolution[s], master_->subsolution_[sindex]);
			subprimobj += workerlb->subprobs_[s]->getPrimalObjective();
			subdualobj += workerlb->subprobs_[s]->getDualObjective();
			DSPdebugMessage("Scenario %d, primobj %+e, dualobj %+e, primobj_sum %+e, dualobj_sum %+e\n",
					workerlb->subprobs_[s]->sind_,
					workerlb->subprobs_[s]->getPrimalObjective(),
					workerlb->subprobs_[s]->getDualObjective(),
					subprimobj, subdualobj);
		}

		if (parEvalUb_ >= 0 || parFeasCuts_ >= 0 || parOptCuts_ >= 0)
		{
			/** store coupling solutions */
			DSP_RTN_CHECK_THROW(storeCouplingSolutions(coupling_solutions));
			DSPdebugMessage("Stored coupling solutions\n");
			/** generate cuts */
			if (parFeasCuts_ >= 0 || parOptCuts_ >= 0)
			{
#ifdef DSP_HAS_SCIP
				if (model_->isDro()) {
					message_->print(0, "Cut procedures are not implemented for DRO.\n");
				} else {
					DSPdebugMessage("generate Benders cuts\n");
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
					DSPdebugMessage("Benders cuts %d\n", cutsToAdd_->sizeCuts());
					DSPdebug2(cutsToAdd_->printCuts());
				}
#endif
			}
			/** evaluate coupling solutions */
			if (parEvalUb_ >= 0)
			{
				DSPdebugMessage("evaluate coupling solutions\n");
				int bestprimsol = -1;
				double oldub = master_->bestprimobj_;
				for (unsigned i = 0; i < coupling_solutions.size(); ++i)
				{
					if (remainingTime() < 0) break;
					/** set time limit */
					workerub->setTimeLimit(remainingTime());
					/** evaluate upper bounds */
					double newub = workerub->evaluate(coupling_solutions[i]);
					DSPdebugMessage("Current upper bound %+e, time remained %e\n", newub, remainingTime());
					if (newub < master_->bestprimobj_)
					{
						master_->bestprimobj_ = newub;
						bestprimsol = i;
					}
				}
				if (oldub > master_->bestprimobj_)
				{
					itercode_ = 'P';
					DSPdebugMessage("model_->getNumCouplingCols() [%d] <= master_->bestprimsol_.size() [%d]", 
						model_->getNumCouplingCols(), (int)master_->bestprimsol_.size());
					assert(model_->getNumCouplingCols()<=master_->bestprimsol_.size());
					for (int j = 0; j < model_->getNumCouplingCols(); ++j) {
						master_->bestprimsol_[j] = (*coupling_solutions[bestprimsol])[j];
					}
				}
			}
			/** clear stored solutions */
			coupling_solutions.clear();
		}

		/** update problem */
#ifdef DSP_DEBUG		
		DSPdebugMessage("subdualobj_:\n");
		DspMessage::printArray(master_->subdualobj_.size(), &(master_->subdualobj_[0]));
#endif
		double olddual = master_->bestdualobj_;
		master_->updateProblem();
		if (olddual < master_->bestdualobj_) {
			itercode_ = itercode_ == 'P' ? 'B' : 'D';
		}

		/** STOP with iteration limit */
		DSPdebugMessage("Check iteration limit\n");
		if (itercnt_ >= master_->getParPtr()->getIntParam("DD/ITER_LIM"))
		{
			message_->print(1, "The iteration limit is reached.\n");
			master_->status_ = DSP_STAT_LIM_ITERorTIME;
			break;
		}

		/** STOP with time limit */
		DSPdebugMessage("Check time limit\n");
		if (remainingTime() < 1.0)
		{
			message_->print(1, "The time limit (%.2f) is reached.\n", parTimeLimit_);
			master_->status_ = DSP_STAT_LIM_ITERorTIME;
			break;
		}

		/** solve problem */
		double tic = CoinGetTimeOfDay();
		DSP_RTN_CHECK_THROW(master_->solve());
		DSPdebugMessage("Solved the master (%.2f sec).\n", CoinGetTimeOfDay() - tic);

		printIterInfo();

		/** termination test */
		if (master_->terminationTest() == DSP_STAT_MW_STOP)
		{
			master_->status_ = DSP_STAT_OPTIMAL;
			break;
		}
		/** check duality gap */
		if (master_->getRelDualityGap() < master_->getParPtr()->getDblParam("DD/STOP_TOL"))
		{
			message_->print(1, "The duality gap limit %+e is reached.\n", master_->getRelDualityGap());
			master_->status_ = DSP_STAT_STOPPED_GAP;
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
		if (model_->isStochastic()) {
			if (model_->isDro()) {
				Ps = master_primsol + master_->getNumCols() - tss->getNumScenarios();
			} else {
				Ps = tss->getProbability();
			}
		}
		master_primsol = NULL;

		/** update subproblems */
		int sindex = -1;
		double probability = -1.0;
		for (int s = 0; s < model_->getNumSubproblems(); ++s)
		{
			sindex = workerlb->subprobs_[s]->sind_;
			workerlb->subprobs_[s]->theta_ = thetas[sindex];
			if (model_->isStochastic()) {
				probability = Ps[sindex];
			}
			// DSPdebugMessage("s = %d, probability = %e, lambdas = \n", s, probability);
			// DspMessage::printArray(model_->getNumSubproblemCouplingRows(sindex), lambdas[sindex]);
			workerlb->subprobs_[s]->updateProblem(lambdas[sindex], probability, master_->bestprimobj_);
			/** apply Benders cuts */
			if (parFeasCuts_ >= 0 || parOptCuts_ >= 0)
			{
				workerlb->subprobs_[s]->pushCuts(cutsToAdd_);
				DSPdebugMessage("Pushed %d Benders cuts.\n", cutsToAdd_->sizeCuts());
			}
		}
	}

	if (parEvalUb_ >= 0 && model_->isStochastic()) {
		DdWorkerUB * workerub = NULL;
		for (unsigned i = 0; i < worker_.size(); ++i)
			if (worker_[i]->getType() == DdWorker::UB) {
				workerub = dynamic_cast<DdWorkerUB*>(worker_[i]);
				break;
			}
		
		/** evaluate UB to get primal solution for each scenario */
		workerub->evaluate(model_->getNumCouplingCols(), &master_->bestprimsol_[0]);
		for (int s = 0; s < tss->getNumScenarios(); ++s) {
			CoinCopyN(workerub->primsols_[s].data(), tss->getNumCols(1), 
				&master_->bestprimsol_[tss->getNumCols(0) + s * tss->getNumCols(1)]);
			//DspMessage::printArray(tss->getNumCols(1), workerub->primsols_[s]);
		}
		DSPdebugMessage2("primsol_:\n");
		DSPdebug2(DspMessage::printArray(model_->getFullModelNumCols(), master_->bestprimsol_.data()));
	}

	/** release shallow-copy of pointers */
	for (int i = 0; i < model_->getNumSubproblems(); ++i)
		lambdas[i] = NULL;

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

#ifdef DSP_HAS_SCIP

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
		DSPdebug(message_->printArray(solutions[i]));

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
					DSPdebug2(rc->print());
					cuts.insertIfNotDuplicate(*rc);

					// TODO: this causes an infinite loop.
					//ret = DSP_STAT_MW_RESOLVE;
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

#endif
