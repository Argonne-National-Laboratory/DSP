/*
 * DdMWSync.cpp
 *
 *  Created on: Mar 23, 2016
 *      Author: kibaekkim
 */

//#define DSP_DEBUG
#include "Model/TssModel.h"
#include "Solver/DualDecomp/DdMWSync.h"
#include "Solver/DualDecomp/DdMasterTr.h"
#ifdef DSP_HAS_OOQP
#include "Solver/DualDecomp/DdMasterDsb.h"
#endif
#include "Solver/DualDecomp/DdMasterSubgrad.h"

DdMWSync::DdMWSync(
		MPI_Comm     comm,   /**< MPI communicator */
		DecModel *   model,  /**< model pointer */
		DspParams *  par,    /**< parameters */
		DspMessage * message /**< message pointer */):
DdMWPara(comm,model,par,message) {}

DdMWSync::DdMWSync(const DdMWSync& rhs) :
DdMWPara(rhs) {}

DdMWSync::~DdMWSync() {}

DSP_RTN_CODE DdMWSync::init()
{
	BGN_TRY_CATCH

	/** This should be before DdMwPara::init(); */
	sync_ = true;

	/** initialize MPI communication settings */
	DdMWPara::init();

	if (comm_rank_ == 0)
	{
		/** create master */
		switch (par_->getIntParam("DD/MASTER_ALGO"))
		{
		case Simplex:
		case IPM:
		case IPM_Feasible:
			master_ = new DdMasterTr(model_, par_, message_);
			break;
// #ifdef DSP_HAS_OOQP
// 		case DSBM:
// 			master_ = new DdMasterDsb(model_, par_, message_);
// 			break;
// #endif
		case Subgradient:
			master_ = new DdMasterSubgrad(model_, par_, message_);
			break;
		}
		/** initialize master */
		master_->init();
	}
	else
	{
		/** clear time log */
		time_lb_.clear();

		/** create workers */
		if (lb_comm_rank_ >= 0)
		{
			/** create LB worker */
			DSPdebugMessage("Rank %d creates a worker for lower bounds.\n", comm_rank_);
			worker_.push_back(new DdWorkerLB(model_, par_, message_));
			/** create CG worker */
			if (parFeasCuts_ >= 0 || parOptCuts_ >= 0)
			{
#ifdef DSP_HAS_SCIP
				DSPdebugMessage("Rank %d creates a worker for Benders cut generation.\n", comm_rank_);
				worker_.push_back(new DdWorkerCGBd(model_, par_, message_));
#endif
			}
			/** create UB worker */
			if (parEvalUb_ >= 0)
			{
				DSPdebugMessage("Rank %d creates a worker for upper bounds.\n", comm_rank_);
				worker_.push_back(new DdWorkerUB(model_, par_, message_));
			}
		}
		/** initialize workers */
		for (unsigned i = 0; i < worker_.size(); ++i)
			worker_[i]->init();
	}

        /** reset iteration info */
        itercnt_   = 0;
        iterstime_ = CoinGetTimeOfDay();

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMWSync::finalize()
{
	BGN_TRY_CATCH

#if 0
	char filename[64];
	const char * output_prefix = par_->getStrParam("OUTPUT/PREFIX").c_str();

	/** free master */
	if (master_)
	{
		sprintf(filename, "%s%d-DD.out", output_prefix, comm_rank_);
		writeIterInfo(filename);
		sprintf(filename, "%s%d-Master.out", output_prefix, comm_rank_);
		master_->write(filename);
		master_->finalize();
		FREE_PTR(master_);
	}

	/** free workers */
	for (unsigned i = 0; i < worker_.size(); ++i)
	{
		switch (worker_[i]->getType())
		{
		case DdWorker::LB:
			sprintf(filename, "%s%d-LB.out", output_prefix, comm_rank_);
			break;
		case DdWorker::UB:
			sprintf(filename, "%s%d-UB.out", output_prefix, comm_rank_);
			break;
		case DdWorker::CGBd:
			sprintf(filename, "%s%d-CG.out", output_prefix, comm_rank_);
			break;
		}
		worker_[i]->write(filename);
		worker_[i]->finalize();
		FREE_PTR(worker_[i]);
	}
	worker_.clear();
#else

	/** print lower bounding times */
	if (par_->getBoolParam("DD/LOG_LB_TIME")) {
		if (lb_comm_rank_ >= 0) {
			if (lb_comm_rank_ == 0) {
				printf("\n## Lower bounding time ##\n");
				printf("%f", time_lb_[0]);
				for (unsigned j = 1; j < time_lb_.size(); ++j)
					printf(",%f", time_lb_[j]);
				printf("\n");
			}

			double* time_lb = new double [time_lb_.size()];
			for (int i = 1; i < lb_comm_size_; ++i) {
				if (lb_comm_rank_ == 0) {
					MPI_Recv(time_lb, (int) time_lb_.size(), MPI_DOUBLE, MPI_ANY_SOURCE, 9999, lb_comm_, MPI_STATUS_IGNORE);
					printf("%f", time_lb[0]);
					for (unsigned j = 1; j < time_lb_.size(); ++j)
						printf(",%f", time_lb[j]);
					printf("\n");
				} else if (lb_comm_rank_ == i)
					MPI_Send(&time_lb_[0], (int) time_lb_.size(), MPI_DOUBLE, 0, 9999, lb_comm_);

				MPI_Barrier(lb_comm_);
			}
			delete [] time_lb;
			time_lb = NULL;

			if (lb_comm_rank_ == 0)
				printf("## End of lower bounding time ##\n\n");
		}
		MPI_Barrier(comm_);
	}

	/** free master */
	if (master_)
	{
		master_->finalize();
		FREE_PTR(master_);
	}

	/** free workers */
	for (unsigned i = 0; i < worker_.size(); ++i)
	{
		worker_[i]->finalize();
		FREE_PTR(worker_[i]);
	}
	worker_.clear();
#endif

	/** finalize MPI settings */
	DdMWPara::finalize();

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMWSync::runMaster()
{
	if (comm_rank_ != 0)
		return DSP_RTN_OK;

#define FREE_MEMORY \
	FREE_ARRAY_PTR(sendbuf) \
	FREE_ARRAY_PTR(scounts) \
	FREE_ARRAY_PTR(sdispls) \
	FREE_ARRAY_PTR(recvbuf) \
	FREE_ARRAY_PTR(rcounts) \
	FREE_ARRAY_PTR(rdispls) \
	FREE_ARRAY_PTR(lambdas) \
	thetas = NULL;

	/** MPI_Scatterv message:
	 *   [for each subproblem]
	 *   1 theta
	 *   2 lambda
	 */
	int size_of_sendbuf;     /**< MPI_Scatterv: size of send buffer pointer */
	double * sendbuf = NULL; /**< MPI_Scatterv: send buffer */
	int *    scounts = NULL; /**< MPI_Scatterv: send buffer size for each process */
	int *    sdispls = NULL; /**< MPI_Scatterv: send buffer displacement for each process*/

	/** MPI_Gatherv message:
	 *   [for each subproblem]
	 *   1 subproblem index
	 *   2 primal objective
	 *   3 dual objective
	 *   4 coupling column part of the solution
	 */
	int size_of_recvbuf;     /**< MPI_Gatherv: size of receive buffer pointer */
	double * recvbuf = NULL; /**< MPI_Gatherv: receive buffer */
	int *    rcounts = NULL; /**< MPI_Gatherv: receive buffer size for each process */
	int *    rdispls = NULL; /**< MPI_Gatherv: receive buffer displacement for each process*/

	const double * thetas  = NULL; /**< of master problem */
	double **      lambdas = NULL; /**< of master problem */

	Solutions stored; /**< coupling solutions */
	Solutions dummy_solutions;
	std::vector<double> dummy_double_array;

	/** Benders cuts */
	int ncuts = 0;
	OsiCuts cuts, emptycuts;
	int cg_status = DSP_STAT_MW_CONTINUE;

	MPI_Status status;

	BGN_TRY_CATCH

	/** collect timing results */
	double mt_total = 0.0; /**< total time */
	double mt_solve = 0.0; /**< solution time */
	double mt_idle = 0.0; /**< idle time (due to subproblem solutions) */
	double mts_total = CoinGetTimeOfDay();
	double mts_solve, mts_idle;

	int signal = DSP_STAT_MW_CONTINUE;   /**< signal to communicate with workers */
	DdMaster * master = dynamic_cast<DdMaster*>(master_);

	/** allocate memory for buffer sizes and displacements */
	scounts = new int [subcomm_size_];
	sdispls = new int [subcomm_size_];
	rcounts = new int [subcomm_size_];
	rdispls = new int [subcomm_size_];

	/** initialize send buffer size and displacement,
	 * and calculate size of send buffer pointer */
	size_of_sendbuf = 0;
	for (int i = 0; i < subcomm_size_; ++i)
	{
		scounts[i] = 0;
		for (int j = 0; j < nsubprobs_[i]; ++j)
			scounts[i] += 1 + model_->getNumSubproblemCouplingRows(subprob_indices_[subprob_displs_[i]+j]);
		sdispls[i] = i == 0 ? 0 : sdispls[i-1] + scounts[i-1];
		size_of_sendbuf += scounts[i];
	}

	/** initialize receive buffer size and displacement,
	 * and calculate size of receive buffer pointer */
	size_of_recvbuf = 0;
	rcounts[0] = 0;
	rdispls[0] = 0;
	for (int i = 1; i < subcomm_size_; ++i)
	{
		rcounts[i] = 0;
		for (int j = 0; j < nsubprobs_[i]; ++j)
			rcounts[i] += 3 + model_->getNumSubproblemCouplingCols(i);
		rdispls[i] = rdispls[i-1] + rcounts[i-1];
		size_of_recvbuf += rcounts[i];
	}

	/** allocate memory for message buffers */
	sendbuf = new double [size_of_sendbuf];
	recvbuf = new double [size_of_recvbuf];

	/** allocate memory for lambdas */
	lambdas = new double * [model_->getNumSubproblems()];

	printHeaderInfo();

	/**
	 * This is the main loop to iteratively solve master problem.
	 */
	while (1)
	{
		/** reset cut generation status */
		cg_status = DSP_STAT_MW_CONTINUE;

		/** reset iteration code */
		itercode_ = ' ';

		/** receive message */
		mts_idle = CoinGetTimeOfDay();
		DSPdebugMessage("Rank %d calls MPI_Gatherv.\n", comm_rank_);
		MPI_Gatherv(NULL, 0, MPI_DOUBLE, recvbuf, rcounts, rdispls, MPI_DOUBLE, 0, subcomm_);
		mt_idle += CoinGetTimeOfDay() - mts_idle;

		DSPdebugMessage2("master receive buffer:\n");
		DSPdebug2(for (int i = 0; i < subcomm_size_; ++i) {
			DSPdebugMessage("  rank %d:\n", i);
			message_->printArray(rcounts[i], recvbuf + rdispls[i]);
		});

		/** apply receive message */
		for (int i = 0, pos = 0; i < subcomm_size_; ++i)
		{
			DSPdebugMessage("message count for rank %d: %d\n", i, rcounts[i]);
			for (int s = 0; s < nsubprobs_[i]; ++s)
			{
				int sindex = static_cast<int>(recvbuf[pos++]);
				master->subprimobj_[sindex] = recvbuf[pos++];
				master->subdualobj_[sindex] = recvbuf[pos++];
				CoinCopyN(recvbuf + pos,
						model_->getNumSubproblemCouplingCols(sindex), master->subsolution_[sindex]);
				pos += model_->getNumSubproblemCouplingCols(sindex);
				DSPdebugMessage("-> master, subprob %d primobj %+e\n", sindex, master->subprimobj_[sindex]);
			}
		}
		//master->worker_ = subcomm_size_ - 1;

		if (parEvalUb_ >= 0 || parFeasCuts_ >= 0 || parOptCuts_ >= 0)
		{
			/** set coupling solutions */
			DSP_RTN_CHECK_THROW(setCouplingSolutions(stored));
			/** clear stored solutions */
			stored.clear();
			/** sync Benders cut information */
			cg_status = syncBendersInfo(stored, cuts);
			/** resolve subproblems? */
			if (cg_status == DSP_STAT_MW_RESOLVE)
			{
				printIterInfo();
				DSPdebugMessage("Rank %d: resolve subproblems.\n", comm_rank_);
				continue;
			}
			/** collect upper bounds */
			if (cg_status == DSP_STAT_MW_CONTINUE)
			  syncUpperbound(dummy_solutions, dummy_double_array);
		}

		/** update problem */
		double olddual = master_->bestdualobj_;
		master_->updateProblem();
		if (olddual < master_->bestdualobj_)
		{
			itercode_ = itercode_ == 'P' ? 'B' : 'D';
			if (master_->getLambda())
				CoinCopyN(master_->getLambda(), model_->getNumCouplingRows(), &master_->bestdualsol_[0]);
		}

		/** STOP with iteration limit */
		if (itercnt_ > master_->getParPtr()->getIntParam("DD/ITER_LIM"))
		{
			signal = DSP_STAT_MW_STOP;
			message_->print(1, "The iteration limit is reached.\n");
			master_->status_ = DSP_STAT_LIM_ITERorTIME;
		}

		/** STOP with time limit */
		if (remainingTime() < 1.0)
		{
			signal = DSP_STAT_MW_STOP;
			message_->print(1, "The time limit (%.2f) is reached.\n", parTimeLimit_);
			master_->status_ = DSP_STAT_LIM_ITERorTIME;
		}

		if (signal != DSP_STAT_MW_STOP)
		{
			/** solve problem */
			double tic = CoinGetTimeOfDay();
			master_->solve();
			mt_solve += CoinGetTimeOfDay() - tic;
			DSPdebugMessage("Rank %d solved the master (%.2f sec).\n", comm_rank_, CoinGetTimeOfDay() - tic);

			/** termination test */
			signal = master_->terminationTest();
			if (signal == DSP_STAT_MW_STOP)
			{
				master_->status_ = DSP_STAT_OPTIMAL;
			}
			else if (master_->getRelDualityGap() < master_->getParPtr()->getDblParam("DD/STOP_TOL"))
			{
				signal = DSP_STAT_MW_STOP;
				message_->print(1, "The duality gap limit %+e is reached.\n", master_->getRelDualityGap());
				master_->status_ = DSP_STAT_STOPPED_GAP;
			}
		}
		printIterInfo();

		/** broadcast signal */
		MPI_Bcast(&signal, 1, MPI_INT, 0, comm_);
		if (signal == DSP_STAT_MW_STOP) break;

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

		/** create send buffer */
		for (int i = 0, pos = 0; i < subcomm_size_; ++i)
		{
			for (int j = 0; j < nsubprobs_[i]; ++j)
			{
				int subprob_index = subprob_indices_[subprob_displs_[i]+j];
				sendbuf[pos++] = thetas[subprob_index];
				CoinCopyN(lambdas[subprob_index],
						model_->getNumSubproblemCouplingRows(subprob_index), sendbuf + pos);
				pos += model_->getNumSubproblemCouplingRows(subprob_index);
			}
		}

		DSPdebugMessage2("master send buffer:\n");
		DSPdebug2(for (int i = 0; i < subcomm_size_; ++i) {
			DSPdebugMessage2("  rank %d:\n", i);
			message_->printArray(scounts[i], sendbuf + sdispls[i]);
		});

		/** scatter message */
		if (parEvalUb_ >= 0)
			MPI_Bcast(&(master_->bestprimobj_), 1, MPI_DOUBLE, 0, subcomm_);
		MPI_Scatterv(sendbuf, scounts, sdispls, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, 0, subcomm_);
	}

	if (parEvalUb_ >= 0 && model_->isStochastic()) {
		TssModel* tss = dynamic_cast<TssModel*>(model_);
		double *bestcouplingsol = NULL;

		/** broadcast best primal coupling solution */
		MPI_Bcast(&master_->bestprimsol_[0], model_->getNumCouplingCols(), MPI_DOUBLE, 0, comm_);

		/** gather primal solution for each scenario */
		bestcouplingsol = new double [tss->getNumCols(1) * tss->getNumScenarios()];
		for (int i = 0; i < subcomm_size_; ++i) {
			rcounts[i] = tss->getNumCols(1) * nsubprobs_[i];
			if (i == 0)
				rdispls[i] = 0;
			else
				rdispls[i] = rdispls[i-1] + rcounts[i-1];
		}
		MPI_Gatherv(NULL, 0, MPI_DOUBLE, bestcouplingsol, rcounts, rdispls, MPI_DOUBLE, 0, comm_);

		/** rearrange primal solution */
		for (int i = 0, j = 0; i < subcomm_size_; ++i)
			for (int s = 0; s < nsubprobs_[i]; ++s) {
				CoinCopyN(bestcouplingsol + j * tss->getNumCols(1), tss->getNumCols(1), 
					&master_->bestprimsol_[tss->getNumCols(0) + subprob_indices_[j] * tss->getNumCols(1)]);
				j++;
			}

		FREE_ARRAY_PTR(bestcouplingsol);

		DSPdebugMessage2("primsol_:\n");
		DSPdebug2(DspMessage::printArray(model_->getFullModelNumCols(), &master_->bestprimsol_[0]));
	}

	/** release shallow-copy of pointers */
	for (int i = 0; i < model_->getNumSubproblems(); ++i)
		lambdas[i] = NULL;

	/** get total time spend in the master */
	mt_total += CoinGetTimeOfDay() - mts_total;
	message_->print(0, "Master timing results: total %.2f, solve %.2f, idle %.2f\n", mt_total, mt_solve, mt_idle);

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

DSP_RTN_CODE DdMWSync::runWorker()
{
	if (comm_rank_ == 0)
		return DSP_RTN_OK;

#define FREE_MEMORY         \
	FREE_ARRAY_PTR(sendbuf) \
	FREE_ARRAY_PTR(recvbuf)

#define SIG_BREAK if (signal == DSP_STAT_MW_STOP) break;

	/** MPI_Gatherv message:
	 *   [for each subproblem]
	 *   1 subproblem index
	 *   2 primal objective
	 *   3 dual objective
	 *   4 coupling column part of the solution
	 */
	double * sendbuf = NULL; /**< MPI_Gatherv: send buffer */
	int      scount  = 0;    /**< MPI_Gatherv: send buffer size */

	/** MPI_Scatterv message:
	 *   [for each subproblem]
	 *   1 theta
	 *   2 lambda
	 */
	double * recvbuf = NULL; /**< MPI_Scatterv: receive buffer */
	int      rcount  = 0;    /**< MPI_Scatterv: receive buffer size */

	OsiCuts cuts, emptycuts;
	bool resolveSubprob = false;
	DSP_RTN_CODE cg_status = DSP_STAT_MW_CONTINUE;

	BGN_TRY_CATCH

	/** timing results */
	double st_total = 0.0;
	double st_cg = 0.0;
	double st_ub = 0.0;
	double st_idle = 0.0;
	double sts_total = CoinGetTimeOfDay();
	double sts_lb, sts_cg, sts_ub, sts_idle;

	int signal = DSP_STAT_MW_CONTINUE;
	int narrprocidx  = par_->getIntPtrParamSize("ARR_PROC_IDX"); /**< number of subproblems */
	int * arrprocidx = par_->getIntPtrParam("ARR_PROC_IDX");     /**< subproblem indices */

	/** calculate size of send buffer */
	for (int i = 0; i < narrprocidx; ++i)
		scount += 3 + model_->getNumSubproblemCouplingCols(arrprocidx[i]);

	/** calculate size of receive buffer */
	for (int i = 0; i < narrprocidx; ++i)
		rcount += 1 + model_->getNumSubproblemCouplingRows(arrprocidx[i]);

	/** allocate memory for message buffers */
	sendbuf = new double [scount];
	recvbuf = new double [rcount];

	/** retrieve DdWorkerLB */
	DdWorkerLB * workerlb = NULL;
	if (lb_comm_ != MPI_COMM_NULL)
	{
		assert(worker_[0]->getType()==DdWorker::LB);
		workerlb = dynamic_cast<DdWorkerLB*>(worker_[0]);
		DSPdebugMessage("Rank %d runs DdWorkerLB.\n", comm_rank_);
	}

	/** solutions to derive Benders cuts and evaluate upper bounds */
	Solutions solutions;

	/** loop until when the master signals stop */
	while (1)
	{
		if (lb_comm_ != MPI_COMM_NULL)
		{
			/** set time limit */
			workerlb->setTimeLimit(remainingTime());

			/** Solve subproblems assigned to each process  */
			sts_lb = CoinGetTimeOfDay();
			workerlb->solve();
			time_lb_.push_back(CoinGetTimeOfDay() - sts_lb);

			/** create send buffer */
			for (int s = 0, pos = 0; s < narrprocidx; ++s)
			{
				sendbuf[pos++] = static_cast<double>(workerlb->subprobs_[s]->sind_);
				sendbuf[pos++] = workerlb->subprobs_[s]->getPrimalObjective();
				sendbuf[pos++] = workerlb->subprobs_[s]->getDualObjective();
				CoinCopyN(workerlb->subprobs_[s]->getSiPtr()->getColSolution(), workerlb->subprobs_[s]->ncols_coupling_, sendbuf + pos);
				pos += model_->getNumSubproblemCouplingCols(workerlb->subprobs_[s]->sind_);
				DSPdebugMessage("MW -> worker %d, subprob %d primobj %+e dualobj %+e\n",
						comm_rank_, workerlb->subprobs_[s]->sind_,
						workerlb->subprobs_[s]->getPrimalObjective(), workerlb->subprobs_[s]->getDualObjective());
			}
#ifdef DSP_DEBUG2
			for (int i = 0; i < lb_comm_size_; ++i)
			{
				if (i == lb_comm_rank_)
				{
					DSPdebugMessage("Rank %d: Worker send message (%d):\n", comm_rank_, scount);
					DSPdebug(message_->printArray(scount, sendbuf));
				}
				MPI_Barrier(lb_comm_);
			}
#endif
			/** send message to the master */
			MPI_Gatherv(sendbuf, scount, MPI_DOUBLE, NULL, NULL, NULL, MPI_DOUBLE, 0, subcomm_);
		}

		/** We may generate Benders-type cuts and find upper bounds */
		if (parEvalUb_ >= 0 || parFeasCuts_ >= 0 || parOptCuts_ >= 0)
		{
			/** get coupling solutions */
			DSP_RTN_CHECK_THROW(bcastCouplingSolutions(solutions));

			/** sync Benders cut information */
			sts_cg = CoinGetTimeOfDay();
			cg_status = syncBendersInfo(solutions, cuts);
			st_cg += CoinGetTimeOfDay() - sts_cg;

			/** calculate and sync upper bounds */
			if (cg_status == DSP_STAT_MW_CONTINUE)
			{
				sts_ub = CoinGetTimeOfDay();
				vector<double> upperbounds;
				DSP_RTN_CHECK_THROW(calculateUpperbound(solutions, upperbounds));
				DSP_RTN_CHECK_THROW(syncUpperbound(solutions, upperbounds));
				st_ub += CoinGetTimeOfDay() - sts_ub;
			}

			/** free solutions */
			for (unsigned i = 0; i < solutions.size(); ++i)
				FREE_PTR(solutions[i]);
			solutions.clear();
		}

		/** receive signal from the master */
		MPI_Bcast(&signal, 1, MPI_INT, 0, comm_);
		DSPdebugMessage2("Rank %d received signal %d.\n", comm_rank_, signal);
		SIG_BREAK;

		if (lb_comm_ != MPI_COMM_NULL)
		{
			/** move cuts to a global pool */
			for (int i = 0; i < cuts.sizeCuts(); ++i)
			{
				OsiRowCut * rc = cuts.rowCutPtr(i);
				cutsToAdd_->insert(rc);
			}
			DSPdebug(cuts.printCuts());
			cuts.dumpCuts();

			/** update subproblems */
			if (cg_status == DSP_STAT_MW_CONTINUE)
			{
				double bestprimalobj = COIN_DBL_MAX;
				if (parEvalUb_ >= 0)
					/** receive upper bound */
					MPI_Bcast(&bestprimalobj, 1, MPI_DOUBLE, 0, subcomm_);
				/** receive message from the master */
				sts_idle = CoinGetTimeOfDay();
				MPI_Scatterv(NULL, NULL, NULL, MPI_DOUBLE, recvbuf, rcount, MPI_DOUBLE, 0, subcomm_);
				st_idle += CoinGetTimeOfDay() - sts_idle;
				DSPdebugMessage2("Worker received message (%d):\n", rcount);
				DSPdebug2(message_->printArray(rcount, recvbuf));
				/** parse message */
				for (int s = 0, pos = 0; s < narrprocidx; ++s)
				{
					workerlb->subprobs_[s]->theta_ = recvbuf[pos++];
					workerlb->subprobs_[s]->updateProblem(recvbuf + pos, bestprimalobj);
					/** apply Benders cuts */
					if (parFeasCuts_ >= 0 || parOptCuts_ >= 0)
					{
						workerlb->subprobs_[s]->pushCuts(cutsToAdd_);
						DSPdebugMessage("Rank %d pushed %d Benders cuts.\n", comm_rank_, cutsToAdd_->sizeCuts());
					}
					pos += model_->getNumSubproblemCouplingRows(workerlb->subprobs_[s]->sind_);
				}
			}
		}
	}

	if (parEvalUb_ >= 0 && model_->isStochastic()) {
		TssModel* tss = dynamic_cast<TssModel*>(model_);
		DdWorkerUB * workerub = NULL;
		double *bestcouplingsol = NULL;

		/** broadcast best primal coupling solution */
		bestcouplingsol = new double [model_->getNumCouplingCols()];
		MPI_Bcast(bestcouplingsol, model_->getNumCouplingCols(), MPI_DOUBLE, 0, comm_);

		/** get WorkerUB pointer */
		for (unsigned i = 0; i < worker_.size(); ++i)
			if (worker_[i]->getType() == DdWorker::UB) {
				workerub = dynamic_cast<DdWorkerUB*>(worker_[i]);
				break;
			}

		/** evaluate UB to get primal solution for each scenario */
		double ub = workerub->evaluate(model_->getNumCouplingCols(), bestcouplingsol);

		/** send primal solution to root */
		FREE_ARRAY_PTR(sendbuf);
		int scount = tss->getNumCols(1) * par_->getIntPtrParamSize("ARR_PROC_IDX");
		sendbuf = new double [scount];
		for (int s = 0; s < par_->getIntPtrParamSize("ARR_PROC_IDX"); ++s)
			CoinCopyN(&workerub->primsols_[s][0], tss->getNumCols(1), sendbuf + s * tss->getNumCols(1));
		MPI_Gatherv(sendbuf, scount, MPI_DOUBLE, NULL, NULL, NULL, MPI_DOUBLE, 0, comm_);

		FREE_ARRAY_PTR(bestcouplingsol);

		DSPdebugMessage2("primsol_:\n");
		DSPdebug2(DspMessage::printArray(model_->getFullModelNumCols(), &master_->bestprimsol_[0]));
	}

	/** release pointers */
	arrprocidx = NULL;

	st_total += CoinGetTimeOfDay() - sts_total;

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef SIG_BREAK
#undef FREE_MEMORY
}

/** broadcast coupling solutions */
DSP_RTN_CODE DdMWSync::bcastCouplingSolutions(
		Solutions & solutions /**< solutions to broadcast */)
{
#define FREE_MEMORY         \
	FREE_ARRAY_PTR(number_of_elements) \
	FREE_ARRAY_PTR(indices)            \
	FREE_ARRAY_PTR(elements)

	int * number_of_elements = NULL; /**< number of elements per solution */
	int * indices            = NULL; /**< indices of solution vectors */
	double * elements        = NULL; /**< elements of solution vectors */

	MPI_Status status;
	int number_of_solutions; /** number of solutions to receive */

	BGN_TRY_CATCH

	if (comm_rank_ == 0)
		number_of_solutions = solutions.size();
	else
	{
		/** clear the local solution pool */
		for (unsigned i = 0; i < solutions.size(); ++i)
			FREE_PTR(solutions[i]);
		solutions.clear();
	}

	/** broadcast the number of solutions */
	MPI_Bcast(&number_of_solutions, 1, MPI_INT, 0, comm_);

	/** broadcast the number of elements for each solution */
	number_of_elements = new int [number_of_solutions];
	int total_elements = 0;
	if (comm_rank_ == 0)
	{
		for (int i = 0; i < number_of_solutions; ++i)
		{
			number_of_elements[i] = solutions[i]->getNumElements();
			total_elements += number_of_elements[i];
		}
	}
	MPI_Bcast(number_of_elements, number_of_solutions, MPI_INT, 0, comm_);
	if (comm_rank_ > 0)
	{
		for (int i = 0; i < number_of_solutions; ++i)
			total_elements += number_of_elements[i];
		DSPdebugMessage("Rank %d received total_elements %d\n", comm_rank_, total_elements);
	}

	/** receive indices and elements */
	indices = new int [total_elements];
	elements = new double [total_elements];
	if (comm_rank_ == 0)
	{
		for (int i = 0, j = 0; i < number_of_solutions; ++i)
		{
			CoinCopyN(solutions[i]->getIndices(), solutions[i]->getNumElements(), indices + j);
			CoinCopyN(solutions[i]->getElements(), solutions[i]->getNumElements(), elements + j);
			j += solutions[i]->getNumElements();
		}
	}
	MPI_Bcast(indices, total_elements, MPI_INT, 0, comm_);
	MPI_Bcast(elements, total_elements, MPI_DOUBLE, 0, comm_);

	if (comm_rank_ > 0)
	{
		for (int i = 0, j = 0; i < number_of_solutions; ++i)
		{
			solutions.push_back(new CoinPackedVector(number_of_elements[i], indices + j, elements + j));
			j += number_of_elements[i];
		}
		DSPdebugMessage("Rank %d received solutions.size() %lu\n", comm_rank_, solutions.size());
	}

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

DSP_RTN_CODE DdMWSync::calculateUpperbound(
		Solutions& solutions, /**< solutions to evaluate */
		vector<double>&  upperbounds /**< list of upper bounds */)
{
	/** return if there is no solution to evaluate */
	if (solutions.size() == 0) return DSP_RTN_OK;
	if (parEvalUb_ < 0) return DSP_RTN_OK;

	DdWorkerUB * workerub = NULL;

	BGN_TRY_CATCH

	/** retrieve upper bounding worker */
	for (unsigned i = 0; i < worker_.size(); ++i)
	{
		if (worker_[i]->getType() == DdWorker::UB)
		{
			workerub = dynamic_cast<DdWorkerUB*>(worker_[i]);
			DSPdebugMessage("Rank %d works for upper bounds.\n", comm_rank_);
			break;
		}
	}
	if (workerub == NULL)
	{
		DSPdebugMessage("Rank %d does not work for upper bounds.\n", comm_rank_);
		return DSP_RTN_OK;
	}

	/** clear upper bounds */
	upperbounds.reserve(solutions.size());
	upperbounds.clear();

	/** calculate upper bound for each solution */
	for (unsigned i = 0; i < solutions.size(); ++i)
	{
		/** set time limit */
		workerub->setTimeLimit(remainingTime());
		/** evaluate upper bounds */
		double sumprimobj = workerub->evaluate(solutions[i]);
		DSPdebugMessage("Rank %d: solution index %d sumprimobj %e\n", comm_rank_, i, sumprimobj);
		upperbounds.push_back(sumprimobj);
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMWSync::setCouplingSolutions(Solutions& solutions)
{
	BGN_TRY_CATCH

	/** store coupling solutions */
	DSP_RTN_CHECK_THROW(storeCouplingSolutions(solutions));
	/** scatter/send coupling solutions */
	DSP_RTN_CHECK_THROW(bcastCouplingSolutions(solutions));

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMWSync::syncBendersInfo(
		Solutions solutions,
		OsiCuts & cuts)
{
	MPI_Status status;
	int cg_status = DSP_STAT_MW_CONTINUE;
	if (parFeasCuts_ < 0 && parOptCuts_ < 0) return cg_status;

	BGN_TRY_CATCH

	if (comm_rank_ == 0)
	{
		/** receive CG status */
		MPI_Recv(&cg_status, 1, MPI_INT, lb_comm_root_, DSP_MPI_TAG_CGBD, comm_, &status);
	}
	else if (lb_comm_rank_ >= 0)
	{
		/** clear cut pool */
		for (int i = 0; i < cuts.sizeCuts(); ++i)
		{
			OsiRowCut * rc = cuts.rowCutPtr(i);
			FREE_PTR(rc);
		}
		cuts.dumpCuts();
		/** generate Benders cuts */
		cg_status = generateBendersCuts(lb_comm_, lb_comm_rank_, solutions, cuts);
		/** send CG status */
		if (lb_comm_rank_ == 0)
			MPI_Send(&cg_status, 1, MPI_INT, 0, DSP_MPI_TAG_CGBD, comm_);
		/** broadcast cuts */
		MPIbcastOsiCuts(lb_comm_, &cuts);
		for (int i = 0; i < lb_comm_size_; ++i)
		{
			if (i == lb_comm_rank_)
			{
				DSPdebugMessage("Rank %d prints cuts:\n", comm_rank_);
				DSPdebug(cuts.printCuts());
			}
			MPI_Barrier(lb_comm_);
		}
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return cg_status;
}

DSP_RTN_CODE DdMWSync::syncUpperbound(
		Solutions&      solutions,  /**< number of solutions */
		vector<double>& upperbounds /**< list of upper bounds */) {
#define FREE_MEMORY              \
	FREE_ARRAY_PTR(primobjzeros) \
	FREE_ARRAY_PTR(primobjs)     \
	FREE_ARRAY_PTR(sumub)

	if (parEvalUb_ < 0) return DSP_RTN_OK;

	int nsolutions;
	double * primobjzeros = NULL; /**< primal objectives filled with zeros */
	double * primobjs     = NULL; /**< primal objectives from subproblems */
	double * sumub        = NULL; /**< the sum of primal objectives */

	MPI_Status status;

	BGN_TRY_CATCH

	if (comm_rank_ == 0)
	{
		/** receive the number of solutions */
		MPI_Recv(&nsolutions, 1, MPI_INT, lb_comm_root_, DSP_MPI_TAG_UB, comm_, &status);
		DSPdebugMessage("Rank %d nsolutions %d\n", comm_rank_, nsolutions);

		/** allocate memory */
		primobjzeros = new double [comm_size_];
		CoinZeroN(primobjzeros, comm_size_);
	}
	else
	{
		nsolutions = solutions.size();
		/** send the number of solutions */
		if (lb_comm_rank_ == 0) {
			MPI_Send(&nsolutions, 1, MPI_INT, 0, DSP_MPI_TAG_UB, comm_);
		}
	}

	/** allocate memory */
	primobjs = new double [nsolutions];
	
	int bestprimsol = -1;
	if (comm_rank_ == 0)
	{
		double oldprimobj = master_->bestprimobj_;
		/** receive upper bounds */
		MPI_Recv(primobjs, nsolutions, MPI_DOUBLE, lb_comm_root_, DSP_MPI_TAG_UB, comm_, &status);
		/** calculate best primal objective */
		for (int i = 0; i < nsolutions; ++i)
		{
			DSPdebugMessage("solution %d: primal objective %+e\n", i, primobjs[i]);
			if (primobjs[i] < master_->bestprimobj_)
			{
				master_->bestprimobj_ = primobjs[i];
				bestprimsol = i;
			}
		}
		if (oldprimobj > master_->bestprimobj_)
			itercode_ = 'P';

		MPI_Bcast(&bestprimsol, 1, MPI_INT, 0, comm_);

		if (bestprimsol > -1) {
			MPI_Recv(&master_->bestprimsol_[0], model_->getNumCouplingCols(), MPI_DOUBLE, 
				lb_comm_root_, DSP_MPI_TAG_UB, comm_, &status);
			//DspMessage::printArray(model_->getNumCouplingCols(), master_->bestprimsol_);
		}
	}
	else if (lb_comm_ != MPI_COMM_NULL)
	{
		/** take the sum of upper bounds */
		if (lb_comm_rank_ == 0)
			sumub = new double [nsolutions];
		MPI_Reduce(&upperbounds[0], sumub, nsolutions, MPI_DOUBLE, MPI_SUM, 0, lb_comm_);
		/** send the upper bounds to the root */
		if (lb_comm_rank_ == 0) {
			MPI_Send(sumub, nsolutions, MPI_DOUBLE, 0, DSP_MPI_TAG_UB, comm_);
		}
		MPI_Bcast(&bestprimsol, 1, MPI_INT, 0, comm_);
		
		if (bestprimsol > -1 && lb_comm_rank_ == 0) {
			MPI_Send(solutions[bestprimsol]->denseVector(model_->getNumCouplingCols()), model_->getNumCouplingCols(),
				MPI_DOUBLE, 0, DSP_MPI_TAG_UB, comm_);
		}
	}

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMWSync::scatterCouplingSolutions(Solutions & solutions) {

	/** temporary empty solutions */
	Solutions emptySolutions;

	/** return signal */
	int signal = DSP_STAT_MW_CONTINUE;

	BGN_TRY_CATCH

	double t = CoinGetTimeOfDay();
	if (comm_rank_ == 0)
	{
		/** After this call, every processor should have the same solutionsGathered. */
		MPIAllgatherCoinPackedVectors(subcomm_, solutions, emptySolutions);

		/** delete local solutions */
		for (unsigned i = 0; i < emptySolutions.size(); ++i)
			FREE_PTR(emptySolutions[i]);
	}
	else
	{
		/** clear the local solution pool */
		for (unsigned i = 0; i < solutions.size(); ++i)
			FREE_PTR(solutions[i]);
		solutions.clear();

		/** After this call, every processor should have the same solutionsGathered. */
		MPIAllgatherCoinPackedVectors(subcomm_, emptySolutions, solutions);
//		for (unsigned i = 0; i < solutions.size(); ++i)
//		{
//			DSPdebugMessage("solutions[%d]:\n", i);
//			message_->printArray(solutions[i]);
//		}
	}
	DSPdebugMessage2("Rank %d: scattered %lu coupling solutions (%.2f sec).\n", comm_rank_, solutions.size(), CoinGetTimeOfDay() - t);

	END_TRY_CATCH_RTN(;,DSP_STAT_MW_STOP)

	return signal;
}
