/*
 * DdMWSync.cpp
 *
 *  Created on: Mar 23, 2016
 *      Author: kibaekkim
 */

#define DSP_DEBUG
#include "Solver/DualDecomp/DdMWSync.h"
#include "Solver/DualDecomp/DdMasterSync.h"

DdMWSync::DdMWSync(
		MPI_Comm          comm,   /**< MPI communicator */
		DdMaster *        master, /**< master problem */
		vector<DdWorker*> worker  /**< worker for finding lower bounds */):
DdMWPara(comm, master, worker) {}

DdMWSync::~DdMWSync() {}

DSP_RTN_CODE DdMWSync::init()
{
	BGN_TRY_CATCH

	DdMWPara::init();
	sync_ = true;

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

	/** Benders cuts */
	int ncuts = 0;
	OsiCuts cuts, emptycuts;
	int cg_status = DSP_STAT_MW_CONTINUE;

	MPI_Status status;

	BGN_TRY_CATCH

	int signal = DSP_STAT_MW_CONTINUE;   /**< signal to communicate with workers */
	DdMasterSync * master = dynamic_cast<DdMasterSync*>(master_);

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

	/**
	 * 0123456789012345678901234567890123456789
	 * iter      curobj     primobj dualobj gap time
	 *  %4d  %+10e  %+10e  %+10e  %6.2f  %6.1f
	 */
	message_->print(1, "  %4s  %13s  %13s  %13s  %6s  %6s\n",
			"iter", "curobj", "primobj", "dualobj", "gap(%)", "times");

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

		/** receive message */
		DSPdebugMessage("Rank %d calls MPI_Gatherv.\n", comm_rank_);
		MPI_Gatherv(NULL, 0, MPI_DOUBLE, recvbuf, rcounts, rdispls, MPI_DOUBLE, 0, subcomm_);

		DSPdebugMessage2("master receive buffer:\n");
		DSPdebug2(for (int i = 0; i < subcomm_size_; ++i) {
			DSPdebugMessage("  rank %d:\n", i);
			message_->printArray(rcounts[i], recvbuf + rdispls[i]);
		});

		/** apply receive message */
		master->nsubprobs_ = 0;
		for (int i = 0, j = 0, pos = 0; i < subcomm_size_; ++i)
		{
			message_->print(3, "message count for rank %d: %d\n", i, rcounts[i]);
			for (int s = 0; s < nsubprobs_[i]; ++s)
			{
				master->subindex_[j] = static_cast<int>(recvbuf[pos++]);
				master->subprimobj_[j] = recvbuf[pos++];
				master->subdualobj_[j] = recvbuf[pos++];
				CoinCopyN(recvbuf + pos, model_->getNumSubproblemCouplingCols(master->subindex_[j]), master->subsolution_[j]);
				pos += model_->getNumSubproblemCouplingCols(master->subindex_[j]);
				message_->print(5, "-> master, subprob %d primobj %+e\n", master->subindex_[j], master->subprimobj_[j]);
				j++;
			}
			master->nsubprobs_ += nsubprobs_[i];
		}
		master->worker_ = subcomm_size_ - 1;

		/** calculate dual objective */
		double dualobj = 0.0;
		for (int s = 0; s < master->nsubprobs_; ++s)
			dualobj += master->subprimobj_[s];
		if (dualobj > master_->bestdualobj_)
		{
			master_->bestdualobj_ = dualobj;
			itercode_ = '*';
		}

		/** STOP with small gap */
		if (itercnt_ > master_->getParPtr()->getIntParam("ITER_LIM"))
		{
			signal = DSP_STAT_MW_STOP;
			message_->print(1, "STOP with iteration limit.\n");
		}
		else if (master_->getRelDualityGap() < master_->getParPtr()->getDblParam("DD/STOP_TOL"))
		{
			signal = DSP_STAT_MW_STOP;
			message_->print(1, "STOP with gap %+e.\n", master_->getRelDualityGap());
		}

		/** broadcast signal */
		DSPdebugMessage("Rank %d send signal %d.\n", comm_rank_, signal);
		/** TODO This is a really bad synchronization. */
		MPI_Bcast(&signal, 1, MPI_INT, 0, comm_);
		if (signal == DSP_STAT_MW_STOP)
		{
			printIterInfo();
			break;
		}

		if (parEvalUb_ >= 0 || parFeasCuts_ >= 0 || parOptCuts_ >= 0)
		{
			/** set coupling solutions */
			setCouplingSolutions(stored);

			/** sync Benders-cut info */
			cg_status = syncBendersInfo(stored, cuts);
			/** clear stored solutions */
			stored.clear();

			/** collect upper bounds */
			if (cg_status == DSP_STAT_MW_CONTINUE)
				syncUpperbound();

			/** resolve subproblems? */
			if (cg_status == DSP_STAT_MW_RESOLVE)
			{
				printIterInfo();
				DSPdebugMessage("Rank %d: resolve subproblems.\n", comm_rank_);
				continue;
			}
		}

		/** update problem */
		master_->updateProblem();

		/** solve problem */
		double tic = CoinGetTimeOfDay();
		master_->solve();
		DSPdebugMessage("Rank %d solved the master (%.2f sec).\n", comm_rank_, CoinGetTimeOfDay() - tic);

		/** returns continue or stop signal */
		signal = master_->terminationTest();

		printIterInfo();

		/** broadcast signal */
		DSPdebugMessage2("Rank %d send signal %d.\n", comm_rank_, signal);
		/** TODO This is a really bad synchronization. */
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
				CoinCopyN(lambdas[subprob_index], model_->getNumSubproblemCouplingRows(subprob_index), sendbuf + pos);
				pos += model_->getNumSubproblemCouplingRows(subprob_index);
			}
		}

		DSPdebugMessage2("master send buffer:\n");
		DSPdebug2(for (int i = 0; i < subcomm_size_; ++i) {
			DSPdebugMessage2("  rank %d:\n", i);
			message_->printArray(scounts[i], sendbuf + sdispls[i]);
		});

		/** scatter message */
		MPI_Scatterv(sendbuf, scounts, sdispls, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, 0, subcomm_);
	}

	/** set best dual objective */
	master_->bestdualobj_ = master_->primobj_;

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
	DSP_RTN_CODE cg_status;

	BGN_TRY_CATCH

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
	if (comm_color_ == comm_color_main)
	{
		assert(worker_[0]->getType()==DdWorker::LB);
		workerlb = dynamic_cast<DdWorkerLB*>(worker_[0]);
	}

	/** solutions to derive Benders cuts and evaluate upper bounds */
	Solutions solutions;

	/** loop until when the master signals stop */
	while (1)
	{
		if (comm_color_ == comm_color_main)
		{
			/** Solve subproblems assigned to each process  */
			workerlb->solve();
			/** create send buffer */
			for (int s = 0, pos = 0; s < narrprocidx; ++s)
			{
				sendbuf[pos++] = static_cast<double>(workerlb->subprobs_[s]->sind_);
				sendbuf[pos++] = workerlb->subprobs_[s]->getPrimalBound();
				sendbuf[pos++] = workerlb->subprobs_[s]->getDualBound();
				CoinCopyN(workerlb->subprobs_[s]->si_->getSolution(), workerlb->subprobs_[s]->ncols_coupling_, sendbuf + pos);
				pos += model_->getNumSubproblemCouplingCols(workerlb->subprobs_[s]->sind_);
				message_->print(5, "MW -> worker %d, subprob %d primobj %+e dualobj %+e\n",
						comm_rank_, workerlb->subprobs_[s]->sind_,
						workerlb->subprobs_[s]->getPrimalBound(), workerlb->subprobs_[s]->getDualBound());
			}
			DSPdebugMessage2("Rank %d: Worker send message (%d):\n", comm_rank_, scount);
			DSPdebug2(message_->printArray(scount, sendbuf));
			/** send message to the master */
			MPI_Gatherv(sendbuf, scount, MPI_DOUBLE, NULL, NULL, NULL, MPI_DOUBLE, 0, subcomm_);
		}

		/** receive signal from the master */
		MPI_Bcast(&signal, 1, MPI_INT, 0, comm_);
		DSPdebugMessage2("Rank %d received signal %d.\n", comm_rank_, signal);
		SIG_BREAK;

		/** We may generate Benders-type cuts and find upper bounds */
		if (parEvalUb_ >= 0 || parFeasCuts_ >= 0 || parOptCuts_ >= 0)
		{
			/** get coupling solutions */
			signal = getCouplingSolutions(solutions);
			SIG_BREAK;

			/** sync Benders cut information */
			cg_status = syncBendersInfo(solutions, cuts);

			/** calculate and sync upper bounds */
			if (cg_status == DSP_STAT_MW_CONTINUE)
			{
				vector<double> upperbounds;
				calculateUpperbound(solutions, upperbounds);
				syncUpperbound(solutions.size(), upperbounds);
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

		if (comm_color_ == comm_color_main)
		{
			/** move cuts to a global pool */
			for (int i = 0; i < cuts.sizeCuts(); ++i)
			{
				OsiRowCut * rc = cuts.rowCutPtr(i);
				cutsToAdd_->insert(rc);
			}
			cuts.dumpCuts();

			/** apply Benders cuts */
			if (parFeasCuts_ >= 0 || parOptCuts_ >= 0)
			{
				workerlb->subprobs_[0]->pushCuts(cutsToAdd_);
				DSPdebugMessage("Rank %d pushed %d Benders cuts.\n", comm_rank_, cutsToAdd_->sizeCuts());
			}

			/** update subproblems */
			if (cg_status == DSP_STAT_MW_CONTINUE)
			{
				/** receive message from the master */
				MPI_Scatterv(NULL, NULL, NULL, MPI_DOUBLE, recvbuf, rcount, MPI_DOUBLE, 0, subcomm_);
				DSPdebugMessage2("Worker received message (%d):\n", rcount);
				DSPdebug2(message_->printArray(rcount, recvbuf));
				/** parse message */
				for (int s = 0, pos = 0; s < narrprocidx; ++s)
				{
					workerlb->subprobs_[s]->theta_ = recvbuf[pos++];
					workerlb->subprobs_[s]->updateProblem(recvbuf + pos);
					pos += model_->getNumSubproblemCouplingRows(workerlb->subprobs_[s]->sind_);
				}
			}
		}
	}

	/** release pointers */
	arrprocidx = NULL;

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef SIG_BREAK
#undef FREE_MEMORY
}

DSP_RTN_CODE DdMWSync::generateBendersCuts(
		Solutions solutions,
		OsiCuts & cuts) {
#define FREE_MEMORY FREE_ARRAY_PTR(aggcut)

	int ret = DSP_STAT_MW_CONTINUE;

	if (solutions.size() == 0) return ret;
	if (parFeasCuts_ < 0 && parOptCuts_ < 0) return ret;
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

//	/** clean cut pool */
//	for (int i = 0; i < cutsToAdd_->sizeCuts(); ++i)
//	{
//		OsiRowCut * rc = cutsToAdd_->rowCutPtr(i);
//		FREE_PTR(rc);
//	}
//	cutsToAdd_->dumpCuts();

	/** retrieve DdWorkerCGBd */
	DdWorkerCGBd * workercg = NULL;
	for (unsigned i = 0; i < worker_.size(); ++i)
	{
		if (worker_[i]->getType() == DdWorker::CGBd)
		{
			workercg = dynamic_cast<DdWorkerCGBd*>(worker_[i]);
			DSPdebugMessage("Rank %d works for generating Benders cuts.\n", comm_rank_);
			break;
		}
	}

	if (workercg == NULL)
	{
		DSPdebugMessage("Rank %d does not work for generating Benders cuts.\n", comm_rank_);
		return ret;
	}

	/** downcast to TssModel */
	tssmodel = dynamic_cast<TssModel*>(model_);

	/** allocate memory */
	aggcut = new double [tssmodel->getNumCols(0) + 1];

	/** generate cuts */
	for (unsigned i = 0; i < solutions.size(); ++i)
	{
		/** create local cut pools */
		OsiCuts cuts, allcuts;

		/** generate cuts at a solution */
		DSP_RTN_CHECK_THROW(workercg->generateCuts(solutions[i], &cuts, cuttype), "generateCuts", "DdWorkerCGBd");
		DSPdebugMessage("Rank %d: Benders cut generator generated %d cuts for solution %d.\n", comm_rank_, cuts.sizeCuts(), i);

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

		/** MPI_Allreduce */
		bool genFeasibilityCuts = false;
		bool genOptimalityCuts  = false;
		MPI_Allreduce(&hasFeasibilityCuts, &genFeasibilityCuts, 1, MPI_C_BOOL, MPI_LOR, cg_comm_);
		MPI_Allreduce(&hasOptimalityCuts, &genOptimalityCuts, 1, MPI_C_BOOL, MPI_LAND, cg_comm_);
		DSPdebugMessage2("Rank %d: hasFeasibilityCuts = %s, hasOptimalityCuts = %s\n",
				comm_rank_, hasFeasibilityCuts ? "true" : "false", hasOptimalityCuts ? "true" : "false");

		if ((genFeasibilityCuts == true && parFeasCuts_ >= 0) ||
			(genOptimalityCuts == true && parOptCuts_ >= 0))
		{
			if (num_comm_colors_ == 1)
				MPIAllgatherOsiCuts(cg_comm_, cuts, allcuts);
			else
				MPIgatherOsiCuts(cg_comm_, cuts, allcuts);

			if (cg_comm_rank_ == 0 || num_comm_colors_ == 1)
			{
				if (genOptimalityCuts)
				{
					/** construct optimality cut */
					CoinZeroN(aggcut, tssmodel->getNumCols(0) + 1);
					for (int j = 0; j < tssmodel->getNumCols(0); ++j)
						aggcut[j] = -(tssmodel->getObjCore(0)[j]);
					aggrhs = 0.0;

					for (int j = 0; j < allcuts.sizeCuts(); ++j)
					{
						OsiRowCut * rc = allcuts.rowCutPtr(j);
						CoinPackedVector cutrow = rc->row();
						for (int k = 0; k < cutrow.getNumElements(); ++k)
							aggcut[cutrow.getIndices()[k]] += cutrow.getElements()[k];
						aggrhs += rc->lb();
						DSPdebugMessage2("allcuts[%d]:\n", j);
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

				if (genFeasibilityCuts)
				{
					/** copy cut pointers */
					for (int j = 0; j < allcuts.sizeCuts(); ++j)
					{
						OsiRowCut * rc = allcuts.rowCutPtr(j);
						DSPdebug(rc->print());
						cuts.insertIfNotDuplicate(*rc);
						ret = DSP_STAT_MW_RESOLVE;
					}
				}
			}
		}

		/** free cuts */
		for (int j = 0; j < cuts.sizeCuts(); ++j)
		{
			OsiRowCut * rc = cuts.rowCutPtr(j);
			FREE_PTR(rc);
		}
		cuts.dumpCuts();
		for (int j = 0; j < allcuts.sizeCuts(); ++j)
		{
			OsiRowCut * rc = allcuts.rowCutPtr(j);
			FREE_PTR(rc);
		}
		allcuts.dumpCuts();
	}

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return ret;
#undef FREE_MEMORY
}

DSP_RTN_CODE DdMWSync::recvCouplingSolutions(Solutions &solutions) {
#define FREE_MEMORY         \
	FREE_ARRAY_PTR(number_of_elements) \
	FREE_ARRAY_PTR(indices)            \
	FREE_ARRAY_PTR(elements)

	assert(comm_color_ > comm_color_main);
	assert(num_comm_colors_ > 1);

	int * number_of_elements = NULL; /**< number of elements per solution */
	int * indices            = NULL; /**< indices of solution vectors */
	double * elements        = NULL; /**< elements of solution vectors */

	/** return signal */
	int signal = DSP_STAT_MW_CONTINUE;

	MPI_Status status;
	int number_of_solutions; /** number of solutions to receive */

	BGN_TRY_CATCH

	/** clear the local solution pool */
	for (unsigned i = 0; i < solutions.size(); ++i)
		FREE_PTR(solutions[i]);
	solutions.clear();

	DSPdebugMessage("rank %d comm_key %d\n", comm_rank_, comm_key_);

	/** receive the number of solutions */
	for (int k = 1; k < num_comm_colors_; ++k)
		if (comm_root_keys_[k] == comm_rank_)
		{
			MPI_Recv(&number_of_solutions, 1, MPI_INT, 0, DSP_MPI_TAG_SOLS, comm_, &status);
			DSPdebugMessage("Rank %d received the number_of_solutions %d\n", comm_rank_, number_of_solutions);
		}
	if (cg_comm_ != MPI_COMM_NULL)
		MPI_Bcast(&number_of_solutions, 1, MPI_INT, 0, cg_comm_);
	if (ub_comm_ != MPI_COMM_NULL)
		MPI_Bcast(&number_of_solutions, 1, MPI_INT, 0, ub_comm_);

	/** receive the number of elements for each solution */
	number_of_elements = new int [number_of_solutions];
	for (int k = 1; k < num_comm_colors_; ++k)
		if (comm_root_keys_[k] == comm_rank_)
		{
			MPI_Recv(number_of_elements, number_of_solutions, MPI_INT, 0, DSP_MPI_TAG_SOLS, comm_, &status);
		}
	if (cg_comm_ != MPI_COMM_NULL)
		MPI_Bcast(number_of_elements, number_of_solutions, MPI_INT, 0, cg_comm_);
	if (ub_comm_ != MPI_COMM_NULL)
		MPI_Bcast(number_of_elements, number_of_solutions, MPI_INT, 0, ub_comm_);

	int total_elements = 0;
	for (int i = 0; i < number_of_solutions; ++i)
		total_elements += number_of_elements[i];
	DSPdebugMessage("total_elements %d\n", total_elements);

	/** receive indices and elements */
	indices = new int [total_elements];
	elements = new double [total_elements];
	for (int k = 1; k < num_comm_colors_; ++k)
		if (comm_root_keys_[k] == comm_rank_)
		{
			MPI_Recv(indices, total_elements, MPI_INT, 0, DSP_MPI_TAG_SOLS, comm_, &status);
			MPI_Recv(elements, total_elements, MPI_DOUBLE, 0, DSP_MPI_TAG_SOLS, comm_, &status);
			for (int i = 0, pos = 0; i < number_of_solutions; ++i)
			{
				DSPdebugMessage2("solution %d:\n", i);
				DSPdebug2({
					for (int j = 0; j < number_of_elements[i]; ++j)
					{
						if (j > 0 && j % 5 == 0) printf("\n");
						printf("  [%6d] %+e", indices[pos+j], elements[pos+j]);
					}
					printf("\n");
				});
				pos += number_of_elements[i];
			}
		}
	if (cg_comm_ != MPI_COMM_NULL)
	{
		MPI_Bcast(indices, total_elements, MPI_INT, 0, cg_comm_);
		MPI_Bcast(elements, total_elements, MPI_DOUBLE, 0, cg_comm_);
	}
	if (ub_comm_ != MPI_COMM_NULL)
	{
		MPI_Bcast(indices, total_elements, MPI_INT, 0, ub_comm_);
		MPI_Bcast(elements, total_elements, MPI_DOUBLE, 0, ub_comm_);
	}

	for (int i = 0, j = 0; i < number_of_solutions; ++i)
	{
		solutions.push_back(new CoinPackedVector(number_of_elements[i], indices + j, elements + j));
		j += number_of_elements[i];
	}

	DSPdebugMessage("solutions.size() %lu\n", solutions.size());

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_STAT_MW_STOP)

	FREE_MEMORY

	return signal;
#undef FREE_MEMORY
}

DSP_RTN_CODE DdMWSync::sendCouplingSolutions(Solutions solutions)
{
#define FREE_MEMORY \
		FREE_ARRAY_PTR(number_of_elements) \
		FREE_ARRAY_PTR(indices)            \
		FREE_ARRAY_PTR(elements)

	int * number_of_elements = NULL; /**< number of elements per solution */
	int * indices            = NULL; /**< indices of solution vectors */
	double * elements        = NULL; /**< elements of solution vectors */

	BGN_TRY_CATCH

	int signal = DSP_STAT_MW_CONTINUE;

	/** send coupling solutions to root process for cut generation and upper bounding */
	int number_of_solutions = solutions.size();
	for (int k = 1; k < num_comm_colors_; ++k)
	{
		DSPdebugMessage("send message to %d\n", comm_root_keys_[k]);
		MPI_Send(&number_of_solutions, 1, MPI_INT, comm_root_keys_[k], DSP_MPI_TAG_SOLS, comm_);
	}

	/** send the number of elements */
	number_of_elements = new int [number_of_solutions];
	int total_elements = 0;
	for (int i = 0; i < number_of_solutions; ++i)
	{
		number_of_elements[i] = solutions[i]->getNumElements();
		total_elements += number_of_elements[i];
	}
	for (int k = 1; k < num_comm_colors_; ++k)
		MPI_Send(number_of_elements, number_of_solutions, MPI_INT, comm_root_keys_[k], DSP_MPI_TAG_SOLS, comm_);

	/** send indices and elements */
	indices = new int [total_elements];
	elements = new double [total_elements];
	for (int i = 0, j = 0; i < number_of_solutions; ++i)
	{
		CoinCopyN(solutions[i]->getIndices(), solutions[i]->getNumElements(), indices + j);
		CoinCopyN(solutions[i]->getElements(), solutions[i]->getNumElements(), elements + j);
		j += solutions[i]->getNumElements();
	}
	for (int k = 1; k < num_comm_colors_; ++k)
	{
		MPI_Send(indices, total_elements, MPI_INT, comm_root_keys_[k], DSP_MPI_TAG_SOLS, comm_);
		MPI_Send(elements, total_elements, MPI_DOUBLE, comm_root_keys_[k], DSP_MPI_TAG_SOLS, comm_);
	}

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

DSP_RTN_CODE DdMWSync::calculateUpperbound(
		Solutions solutions, /**< solutions to evaluate */
		vector<double>&  upperbounds /**< list of upper bounds */)
{
	/** return if there is no solution to evaluate */
	if (solutions.size() == 0) return DSP_RTN_OK;
	if (parEvalUb_ < 0) return DSP_RTN_OK;
	if (ub_comm_rank_ < 0) return DSP_RTN_OK;

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
		/** fix coupling solutions */
		workerub->fixCouplingVariableValues(solutions[i]);

		/** solve */
		workerub->solve();

		/** take minimum objective */
		double sumprimobj = 0.0;
		for (unsigned s = 0; s < workerub->subprobs_.size(); ++s)
			sumprimobj += workerub->subprobs_[s]->getPrimalBound();
		DSPdebugMessage("Rank %d: solution index %d sumprimobj %e\n", comm_rank_, i, sumprimobj);
		upperbounds.push_back(sumprimobj);
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMWSync::storeCouplingSolutions(Solutions & stored) {

	BGN_TRY_CATCH

	DdMasterSync * master  = dynamic_cast<DdMasterSync*>(master_);

	/** store solutions to distribute */
	for (int s = 0; s < master->nsubprobs_; ++s)
	{
		int nx = model_->getNumSubproblemCouplingCols(master->subindex_[s]);

		DSPdebugMessage2("ubSolutions_ %lu\n", ubSolutions_.size());
		DSPdebugMessage2("solution[%d] nx %d:\n", s, nx);
		DSPdebug2(message_->printArray(nx, master->subsolution_[s]));

		CoinPackedVector * x = duplicateSolution(
				nx, master->subsolution_[s], ubSolutions_);
		if (x != NULL)
		{
			DSPdebugMessage2("Coupling solution:\n");
			DSPdebug2(DspMessage::printArray(nx, master->subsolution_[s]));

			/** store solution */
			ubSolutions_.push_back(x);
			stored.push_back(x);
		}
	}
	DSPdebugMessage("Rank %d stored %lu solutions.\n", comm_rank_, stored.size());

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMWSync::recvBendersCuts(OsiCuts & cuts) {

	if (num_comm_colors_ == 1) return DSP_RTN_OK;
	if (parFeasCuts_ < 0 && parOptCuts_ < 0) return DSP_RTN_OK;

	int recv_flag;
	MPI_Status status;

	BGN_TRY_CATCH

	/** check if there is a message to receive */
	MPI_Iprobe(comm_root_keys_[1], DSP_MPI_TAG_CGBD, comm_, &recv_flag, &status);
	if (recv_flag == 0) return DSP_RTN_OK;

	/** receive cuts */
	MPIrecvOsiCuts(comm_, comm_root_keys_[1], cuts, DSP_MPI_TAG_CGBD);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMWSync::getCouplingSolutions(Solutions& solutions) {
	if (num_comm_colors_ == 1)
		/** scatter coupling solutions */
		return scatterCouplingSolutions(solutions);
	else
		/** receive coupling solutions from the master*/
		return recvCouplingSolutions(solutions);
}

DSP_RTN_CODE DdMWSync::setCouplingSolutions(Solutions& solutions) {
	/** store coupling solutions */
	storeCouplingSolutions(solutions);
	/** scatter/send coupling solutions */
	if (num_comm_colors_ == 1)
		return scatterCouplingSolutions(solutions);
	else
		return sendCouplingSolutions(solutions);
}

DSP_RTN_CODE DdMWSync::syncBendersInfo(
		Solutions solutions,
		OsiCuts & cuts)
{
	MPI_Status status;
	DSP_RTN_CODE cg_status = DSP_STAT_MW_CONTINUE;

	if (solutions.size() <= 0) return cg_status;
	if (parFeasCuts_ < 0 && parOptCuts_ < 0) return cg_status;

	BGN_TRY_CATCH

	switch (num_comm_colors_)
	{
	case 1:
		if (comm_rank_ == 0)
		{
			/** receive cut generation status */
			MPI_Recv(&cg_status, 1, MPI_INT, comm_root_keys_[1], DSP_MPI_TAG_CGBD, comm_, &status);
		}
		else
		{
			/** generate Benders cuts */
			cg_status = generateBendersCuts(solutions, cuts);
			if (cg_comm_rank_ == 0)
			{
				/** send the cut generation status */
				MPI_Send(&cg_status, 1, MPI_INT, 0, DSP_MPI_TAG_CGBD, comm_);
			}
		}
		break;
	case 2:
	case 3:
	{
		/** clear cut pool */
		for (int i = 0; i < cuts.sizeCuts(); ++i)
		{
			OsiRowCut * rc = cuts.rowCutPtr(i);
			FREE_PTR(rc);
		}
		cuts.dumpCuts();

		if (comm_rank_ == 0)
		{
			/** receive cut generation status */
			MPI_Recv(&cg_status, 1, MPI_INT, comm_root_keys_[1], DSP_MPI_TAG_CGBD, comm_, &status);
			/** scatter cut generation status */
			MPI_Scatter(&cg_status, 1, MPI_INT, NULL, 0, MPI_INT, 0, subcomm_);
			/** receive Benders cuts */
			recvBendersCuts(cuts);
			/** scatter Benders cuts to the lower bound workers*/
			MPIscatterOsiCuts(subcomm_, cuts, NULL);
		}
		else
		{
			if (comm_color_ == comm_color_cg)
			{
				/** generate Benders cuts */
				cg_status = generateBendersCuts(solutions, cuts);
				if (cg_comm_rank_ == 0)
				{
					/** send the cut generation status */
					MPI_Send(&cg_status, 1, MPI_INT, 0, DSP_MPI_TAG_CGBD, comm_);
					/** send the cuts to the root */
					MPIsendOsiCuts(comm_, 0, cuts, DSP_MPI_TAG_CGBD);
				}
			}
			else if (comm_color_ == comm_color_main)
			{
				/** scatter cut generation status */
				MPI_Scatter(NULL, 0, MPI_INT, &cg_status, 1, MPI_INT, 0, subcomm_);
				/** receive Benders cuts */
				MPIscatterOsiCuts(subcomm_, OsiCuts(), &cuts);
			}
			/** TODO comm_color_ub? */
		}
		break;
	}
	default:
		throw "Unexpected value for the number of communication colors.";
		break;
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return cg_status;
}

DSP_RTN_CODE DdMWSync::syncUpperbound(
		int nsolutions, /**< number of solutions */
		vector<double> upperbounds /**< list of upper bounds */) {
#define FREE_MEMORY              \
	FREE_ARRAY_PTR(primobjzeros) \
	FREE_ARRAY_PTR(primobjs)     \
	FREE_ARRAY_PTR(sumub)

	if (parEvalUb_ < 0) return DSP_RTN_OK;

	double * primobjzeros = NULL; /**< primal objectives filled with zeros */
	double * primobjs     = NULL; /**< primal objectives from subproblems */
	double * sumub        = NULL; /**< the sum of primal objectives */

	MPI_Status status;

	BGN_TRY_CATCH

	if (comm_rank_ == 0)
	{
		/** receive the number of solutions */
		MPI_Recv(&nsolutions, 1, MPI_INT, comm_root_keys_[2], DSP_MPI_TAG_UB, comm_, &status);
		DSPdebugMessage("Rank %d nsolutions %d\n", comm_rank_, nsolutions);

		/** allocate memory */
		primobjs = new double [nsolutions];
		primobjzeros = new double [comm_size_];
		CoinZeroN(primobjzeros, comm_size_);
	}
	else
	{
		/** number of solutions */
		if (ub_comm_rank_ == 0)
		{
			/** send the number of solutions */
			MPI_Send(&nsolutions, 1, MPI_INT, 0, DSP_MPI_TAG_UB, comm_);
		}
	}

	switch (num_comm_colors_)
	{
	case 1:
		if (comm_rank_ == 0)
		{
			/** reduce the sum of upper bounds */
			MPI_Reduce(primobjzeros, primobjs, nsolutions, MPI_DOUBLE, MPI_SUM, 0, comm_);
		}
		else
		{
			/** reduce the sum of upper bounds */
			MPI_Reduce(&upperbounds[0], NULL, nsolutions, MPI_DOUBLE, MPI_SUM, 0, comm_);
		}
		break;
	case 2:
	case 3:
		if (comm_rank_ == 0)
		{
			/** receive upper bounds */
			MPI_Recv(primobjs, nsolutions, MPI_DOUBLE, comm_root_keys_[2], DSP_MPI_TAG_UB, comm_, &status);
		}
		else
		{
			/** take the sum of upper bounds */
			if (ub_comm_rank_ == 0)
				sumub = new double [nsolutions];
			MPI_Reduce(&upperbounds[0], sumub, nsolutions, MPI_DOUBLE, MPI_SUM, 0, ub_comm_);
			/** send the upper bounds to the root */
			if (ub_comm_rank_ == 0)
				MPI_Send(sumub, nsolutions, MPI_DOUBLE, 0, DSP_MPI_TAG_UB, comm_);
		}
		break;
	default:
		throw "Unexpected value for the number of communication colors.";
		break;
	}

	if (comm_rank_ == 0)
	{
		/** calculate primal objective */
		for (int i = 0; i < nsolutions; ++i)
		{
			DSPdebugMessage("solution %d: primal objective %+e\n", i, primobjs[i]);
			master_->bestprimobj_ = primobjs[i] < master_->bestprimobj_ ? primobjs[i] : master_->bestprimobj_;
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
