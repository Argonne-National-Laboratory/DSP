/*
 * DdMWSync.cpp
 *
 *  Created on: Mar 23, 2016
 *      Author: kibaekkim
 */

//#define DSP_DEBUG
#include "Solver/DualDecomp/DdMWSync.h"
#include "Solver/DualDecomp/DdMasterSync.h"
#include "Solver/DualDecomp/DdMasterTr.h"
#include "Solver/DualDecomp/DdMasterDsb.h"
#include "Solver/DualDecomp/DdMasterSubgrad.h"
#include "Solver/DualDecomp/DdMasterReg.h"

DdMWSync::DdMWSync(
		MPI_Comm     comm,   /**< MPI communicator */
		DecModel *   model,  /**< model pointer */
		DspParams *  par,    /**< parameters */
		DspMessage * message /**< message pointer */):
DdMWPara(comm,model,par,message) {}

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
			master_ = new DdMasterTr(par_, model_, message_, comm_size_);
			break;
		case DSBM:
			master_ = new DdMasterDsb(par_, model_, message_, comm_size_);
			break;
		case Subgradient:
			master_ = new DdMasterSubgrad(par_, model_, message_, comm_size_);
			break;
		case Regularize_Bundle:
			master_ = new DdMasterReg(par_, model_, message_, comm_size_);
			break;
		}
		/** initialize master */
		master_->init();
	}
	else
	{
		/** create workers */
		if (lb_comm_rank_ >= 0)
		{
			/** create LB worker */
			DSPdebugMessage("Rank %d creates a worker for lower bounds.\n", comm_rank_);
			worker_.push_back(new DdWorkerLB(par_, model_, message_));
			/** create CG worker */
			if (parFeasCuts_ >= 0 || parOptCuts_ >= 0)
			{
				DSPdebugMessage("Rank %d creates a worker for Benders cut generation.\n", comm_rank_);
				worker_.push_back(new DdWorkerCGBd(par_, model_, message_));
			}
			/** create UB worker */
			if (parEvalUb_ >= 0)
			{
				DSPdebugMessage("Rank %d creates a worker for upper bounds.\n", comm_rank_);
				worker_.push_back(new DdWorkerUB(par_, model_, message_));
			}
		}
		/** initialize workers */
		for (unsigned i = 0; i < worker_.size(); ++i)
			worker_[i]->init();
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMWSync::finalize()
{
	BGN_TRY_CATCH

	char filename[64];
	const char * output_prefix = par_->getStrParam("OUTPUT/PREFIX").c_str();

	/** free master */
	if (master_)
	{
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

		DSPdebugMessage("master receive buffer:\n");
		DSPdebug(for (int i = 0; i < subcomm_size_; ++i) {
			DSPdebugMessage("  rank %d:\n", i);
			message_->printArray(rcounts[i], recvbuf + rdispls[i]);
		});

		/** apply receive message */
		master->nsubprobs_ = 0;
		for (int i = 0, j = 0, pos = 0; i < subcomm_size_; ++i)
		{
			message_->print(5, "message count for rank %d: %d\n", i, rcounts[i]);
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

		/** STOP with small gap */
		if (itercnt_ > master_->getParPtr()->getIntParam("ITER_LIM"))
		{
			signal = DSP_STAT_MW_STOP;
			message_->print(0, "The iteration limit is reached.\n");
		}
		/** broadcast signal */
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
				syncUpperbound();
		}

		/** update problem */
		double olddual = master_->bestdualobj_;
		master_->updateProblem();
		if (olddual < master_->bestdualobj_)
			itercode_ = itercode_ == 'P' ? 'B' : 'D';

		/** solve problem */
		double tic = CoinGetTimeOfDay();
		master_->solve();
		DSPdebugMessage("Rank %d solved the master (%.2f sec).\n", comm_rank_, CoinGetTimeOfDay() - tic);

		/** termination test */
		signal = master_->terminationTest();
		/** check duality gap */
		if (signal != DSP_STAT_MW_STOP &&
				master_->getRelDualityGap() < master_->getParPtr()->getDblParam("DD/STOP_TOL"))
		{
			signal = DSP_STAT_MW_STOP;
			message_->print(0, "The duality gap limit %+e is reached.\n", master_->getRelDualityGap());
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
		if (parEvalUb_ >= 0)
			MPI_Bcast(&(master_->bestprimobj_), 1, MPI_DOUBLE, 0, subcomm_);
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
	DSP_RTN_CODE cg_status = DSP_STAT_MW_CONTINUE;

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

		/** receive signal from the master */
		MPI_Bcast(&signal, 1, MPI_INT, 0, comm_);
		DSPdebugMessage("Rank %d received signal %d.\n", comm_rank_, signal);
		SIG_BREAK;

		/** We may generate Benders-type cuts and find upper bounds */
		if (parEvalUb_ >= 0 || parFeasCuts_ >= 0 || parOptCuts_ >= 0)
		{
			/** get coupling solutions */
			signal = bcastCouplingSolutions(solutions);
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

		if (lb_comm_ != MPI_COMM_NULL)
		{
			/** move cuts to a global pool */
			for (int i = 0; i < cuts.sizeCuts(); ++i)
			{
				OsiRowCut * rc = cuts.rowCutPtr(i);
				cutsToAdd_->insert(rc);
			}
			cuts.dumpCuts();

			/** update subproblems */
			if (cg_status == DSP_STAT_MW_CONTINUE)
			{
				double bestprimalobj = COIN_DBL_MAX;
				if (parEvalUb_ >= 0)
					/** receive upper bound */
					MPI_Bcast(&bestprimalobj, 1, MPI_DOUBLE, 0, subcomm_);
				/** receive message from the master */
				MPI_Scatterv(NULL, NULL, NULL, MPI_DOUBLE, recvbuf, rcount, MPI_DOUBLE, 0, subcomm_);
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

	/** release pointers */
	arrprocidx = NULL;

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

	/** return signal */
	int signal = DSP_STAT_MW_CONTINUE;

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

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_STAT_MW_STOP)

	FREE_MEMORY

	return signal;
#undef FREE_MEMORY
}

DSP_RTN_CODE DdMWSync::calculateUpperbound(
		Solutions solutions, /**< solutions to evaluate */
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

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMWSync::setCouplingSolutions(Solutions& solutions) {
	/** store coupling solutions */
	storeCouplingSolutions(solutions);
	/** scatter/send coupling solutions */
	return bcastCouplingSolutions(solutions);
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
		MPI_Recv(&nsolutions, 1, MPI_INT, lb_comm_root_, DSP_MPI_TAG_UB, comm_, &status);
		DSPdebugMessage("Rank %d nsolutions %d\n", comm_rank_, nsolutions);

		/** allocate memory */
		primobjzeros = new double [comm_size_];
		CoinZeroN(primobjzeros, comm_size_);
	}
	else
	{
		/** send the number of solutions */
		if (lb_comm_rank_ == 0)
			MPI_Send(&nsolutions, 1, MPI_INT, 0, DSP_MPI_TAG_UB, comm_);
	}

	/** allocate memory */
	primobjs = new double [nsolutions];

	if (comm_rank_ == 0)
	{
		double oldprimobj = master_->bestprimobj_;
		/** receive upper bounds */
		MPI_Recv(primobjs, nsolutions, MPI_DOUBLE, lb_comm_root_, DSP_MPI_TAG_UB, comm_, &status);
		/** calculate best primal objective */
		for (int i = 0; i < nsolutions; ++i)
		{
			DSPdebugMessage("solution %d: primal objective %+e\n", i, primobjs[i]);
			master_->bestprimobj_ = primobjs[i] < master_->bestprimobj_ ? primobjs[i] : master_->bestprimobj_;
		}
		if (oldprimobj > master_->bestprimobj_)
			itercode_ = 'P';
	}
	else if (lb_comm_ != MPI_COMM_NULL)
	{
		/** take the sum of upper bounds */
		if (lb_comm_rank_ == 0)
			sumub = new double [nsolutions];
		MPI_Reduce(&upperbounds[0], sumub, nsolutions, MPI_DOUBLE, MPI_SUM, 0, lb_comm_);
		/** send the upper bounds to the root */
		if (lb_comm_rank_ == 0)
			MPI_Send(sumub, nsolutions, MPI_DOUBLE, 0, DSP_MPI_TAG_UB, comm_);
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
