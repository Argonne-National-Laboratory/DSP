/*
 * DecDdMpi.cpp
 *
 *  Created on: Dec 10, 2014
 *      Author: kibaekkim, ctjandra
 */

//#define DSP_DEBUG
//#define DO_SERIALIZE
#define USE_ROOT_PROCESSOR 1

#include <vector>

/** DSP */
#include "Utility/StoMessage.h"
#include "Utility/StoUtility.h"
#include "Solver/TssBdSub.h"
#include "Solver/DecDdPrimalMaster.h"
#include "Solver/DecDdMasterSubgrad.h"
#include "Solver/DecDdMasterDSB.h"
#include "Solver/DecDdSub.h"
#include "SolverInterface/SolverInterfaceScip.h"
#include "SolverInterface/OoqpEps.h"
#include "Solver/DecDdMpi.h"

/** COIN */
#include "CoinTime.hpp"

using namespace std;

/** default constructor */
DecDdMpi::DecDdMpi(
		MPI_Comm comm,
		string logfile_prefix) :
	DecSolver(),
	comm_(comm),
	num_comm_groups_(-1),
	comm_group_(-1),
	cuts_(NULL),
	numSyncedCuts_(0),
	bdsub_(NULL),
	numSyncedUbSolutions_(0),
	master_(NULL),
	subproblemSpecs_(NULL),
	multipliers_(NULL),
	tic_(MPI::Wtime()),
	iterCnt_(0),
	nInfeasible_(0),
	ctime_start_(0.),
	wtime_start_(0.),
	logfile_prefix_(logfile_prefix),
	parMasterAlgo_(0),
	parFeasCuts_(-1),
	parOptCuts_(-1),
	parEvalUb_(-1),
	parProcIdxSize_(0),
	parProcIdx_(NULL),
	parNumCutsPerIter_(1),
	parLogDualVars_(false),
	parStopTol_(1.0e-5)
{
	MPI_Comm_rank(comm, &comm_rank_);
	MPI_Comm_size(comm, &comm_size_);
	nsubprobsAtRank_ = new int [comm_size_];
}

DecDdMpi::~DecDdMpi()
{
	FREE_PTR(bdsub_);
	FREE_PTR(master_);
	FREE_ARRAY_PTR(nsubprobsAtRank_);
	FREE_ARRAY_PTR(subproblemSpecs_);
	FREE_ARRAY_PTR(multipliers_);
	for (unsigned int i = 0; i < subprobs_.size(); ++i)
		FREE_PTR(subprobs_[i]);
	for (unsigned int i = 0; i < ubSolutions_.size(); ++i)
		FREE_PTR(ubSolutions_[i]);
	for (unsigned int i = 0; i < solutions_.size(); ++i)
		FREE_PTR(solutions_[i]);
	if (cuts_)
	{
		for (int i = 0; i < cuts_->sizeCuts(); ++i)
		{
			OsiRowCut * rc = cuts_->rowCutPtr(i);
			FREE_PTR(rc);
		}
		cuts_->dumpCuts();
	}
	FREE_ARRAY_PTR(parProcIdx_);
}

/** initialize solver */
STO_RTN_CODE DecDdMpi::initialize()
{
	BGN_TRY_CATCH

	int nsubprobs = model_->getNumSubproblems();

	/** read parameters */
	parMasterAlgo_     = par_->getIntParam("DD/MASTER_ALGO");
	parFeasCuts_       = par_->getIntParam("DD/FEAS_CUTS");
	parOptCuts_        = par_->getIntParam("DD/OPT_CUTS");
	parEvalUb_         = par_->getIntParam("DD/EVAL_UB");
	parNumCutsPerIter_ = par_->getIntParam("DD/NUM_CUTS_PER_ITER");
	parLogDualVars_    = par_->getBoolParam("DD/LOG_DUAL_VARS");
	parStopTol_        = par_->getDblParam("DD/STOP_TOL");

	/** cuts and upper bound only supported with nonanticipativity constraints */
	if (!model_->nonanticipativity() || !model_->isStochastic())
	{
		if (parFeasCuts_ >= 0 || parOptCuts_ >= 0)
		{
			printf("Warning: Feasibility and optimality cuts are not supported in a general decomposition model\n");
		}
		if (parEvalUb_ >= 0)
		{
			printf("Warning: Primal heuristic is not supported in a general decomposition model\n");
		}
		parFeasCuts_ = -1;
		parOptCuts_ = -1;
		parEvalUb_ = -1;
	}

	/** if TssDdProcIdxSet_ not set, initialize subproblems per process in a default manner */
#if USE_ROOT_PROCESSOR == 0
	int effective_comm_rank = comm_rank_ - 1;
	int effective_comm_size = comm_size_ - 1;
#else
	int effective_comm_rank = comm_rank_;
	int effective_comm_size = comm_size_;
#endif
	if (parProcIdxSize_ == 0)
	{
		effective_comm_size = CoinMin(effective_comm_size, model_->getNumSubproblems());
		if (effective_comm_rank < effective_comm_size)
		{
			/** round-robin */
			int size = (int) ((model_->getNumSubproblems() - 1) - effective_comm_rank) / effective_comm_size + 1;
			assert(size > 0);
			parProcIdxSize_ = size;
			parProcIdx_ = new int [size];
			int i = 0;
			for (int s = effective_comm_rank; s < model_->getNumSubproblems(); s += effective_comm_size)
				parProcIdx_[i++] = s;
		}
	}

	/** calculate communication groups */
	num_comm_groups_ = (comm_size_ - 1) / nsubprobs + 1;
	if (comm_size_ > nsubprobs
			&& comm_size_ % nsubprobs > 0
			&& comm_size_ % nsubprobs < comm_rank_ % nsubprobs + 1)
		num_comm_groups_ -= 1;
	comm_group_ = comm_rank_ / nsubprobs;
	DSPdebugMessage("-> comm_size %d comm_rank %d num_comm_groups %d comm_group %d\n",
			comm_size_, comm_rank_, num_comm_groups_, comm_group_);

	/** Sanity check:
	 *    The number of master cuts per iteration should not exceed the number of subproblems;
	 *    should not be less than 1. */
	if (parNumCutsPerIter_ > nsubprobs || parNumCutsPerIter_ < 1)
		parNumCutsPerIter_ = nsubprobs;

	/** Root process can print out. */
	if (comm_rank_ > 0) parLogLevel_ = -100;

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
}

/** solve */
STO_RTN_CODE DecDdMpi::solve()
{
	BGN_TRY_CATCH

#if USE_ROOT_PROCESSOR == 0
	if (comm_size_ < 2)
	{
		printf("Error: At least two processors are required.\n");
		status_ = STO_STAT_STOPPED_MPI;
		return STO_RTN_ERR;
	}
#endif

	/** initialize solver */
	initialize();

	/** tic */
	ctime_start_ = CoinCpuTime();
	wtime_start_ = MPI_Wtime();

	/** create master problem */
	createMaster();

	/** create subproblems */
	createSubproblem();

	/** print time elapsed for creating problems */
	if (parLogLevel_ > 0)
		printf("Time elapsed for creating master and subproblems: %.2f seconds\n",
			MPI_Wtime() - wtime_start_);

	while (time_remains_ >= 0.)
	{
		ticToc(); /**< update wall clock */

		/** solve subproblems: may update primal bound. */
		status_ = solveSubproblem();
		if (status_ != STO_RTN_OK)
			break;

		ticToc(); /**< update wall clock */

		/** synchronize master:
		 * 1. send subproblem solutions to master,
		 * 2. update master problem,
		 * 3. and may update dual bound */
		syncMaster();

		/** solve master */
		solveMaster();

		ticToc(); /**< update wall clock */

		/** timestamp per iteration */
		wtime_elapsed_.push_back(MPI::Wtime() - wtime_start_);

		/** terminate? */
		if (doTerminate() || time_remains_ < 0.)
			break;

		/** increase iteration counter */
		iterCnt_++;

		/** synchronize subproblem: scatter Lagrangian multipliers */
		syncSubproblem();

#ifdef DSP_DEBUG
		if (comm_rank_ == 0)
			DSPdebugMessage("Sychronized subproblems (%.2f seconds)\n", MPI::Wtime() - tic_);
#endif
	}

	/** collect results */
	collectResults();

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
}

/** create subproblem */
STO_RTN_CODE DecDdMpi::createSubproblem()
{
	int * displs = NULL;
	vector<int> subproblemsForCut;

	/** SCENARIO DISTRIBUTION:
	 *    scenarios taken care of by the current processor */
	for (int s = 0; s < parProcIdxSize_; ++s)
		subproblemsForCut.push_back(parProcIdx_[s]);

#if USE_ROOT_PROCESSOR == 0
	if (comm_rank_ != 0)
#endif
	{
		if (parFeasCuts_ >= 0 || parOptCuts_ >= 0 || parEvalUb_ >= 0)
		{

			/** these settings are only supported for stochastic models */
			TssModel * tssModel;
			try
			{
				tssModel = dynamic_cast<TssModel *>(model_);
			}
			catch (const std::bad_cast& e)
			{
				printf("Error: Cuts and primal heuristic are not supported in a general decomposition model\n");
				return STO_RTN_ERR;
			}

			/** create Benders subproblem:
			 *   This is used in generating feasibility/optimality cuts
			 *   and evaluating upper bounds. */
			bdsub_ = new TssBdSub(par_);
			bdsub_->scenarios_ = subproblemsForCut;
			bdsub_->loadProblem(tssModel);
		}

		for (int s = 0; s < parProcIdxSize_; ++s)
		{
			int sub = parProcIdx_[s];
			/** create subproblem instance */
			DecDdSub * subprob = new DecDdSub(sub, par_);
			subprobs_.push_back(subprob);

			/** create subproblem */
			subprob->createProblem(model_);
			assert(subprob->si_);

			SolverInterfaceScip * si = dynamic_cast<SolverInterfaceScip*>(subprob->si_);
			if (si)
			 	subprob->si_->setPrintLevel(CoinMin(CoinMax(0, parLogLevel_ - 3), 5));
			else
			 	subprob->si_->setPrintLevel(CoinMax(0, parLogLevel_ - 3));

			/** add cut generator */
			if (parFeasCuts_ >= 0 || parOptCuts_ >= 0)
			{
				/** Should change later */
				subprob->addCutGenerator(bdsub_);
			}
#if 0
			/** add branch rule for known LB */
			subprob->addBranchrule();
#endif
		}
	}

	/** cut pool */
	cuts_ = new OsiCuts;

	/** Let the root know how many scenarios each process use.
	 *  We do not count the number of subproblems for comm_group > 0,
	 *  in which case the subproblem data is used to generate cuts and to evaluate upper bounds. */
	int nsubprobs_proc = comm_group_ > 0 ? 0 : subprobs_.size();
	MPI_Gather(&nsubprobs_proc, 1, MPI::INT, nsubprobsAtRank_, 1, MPI::INT, 0, comm_);

	if (comm_rank_ == 0)
	{
		subproblemSpecs_ = new int [model_->getNumSubproblems()];
		displs = new int [comm_size_];
		displs[0] = 0;
		for (int i = 1; i < comm_size_; ++i)
			displs[i] = displs[i-1] + nsubprobsAtRank_[i-1];
	}

	/** Let the root know which subproblems each process take care of. */
	MPI_Gatherv(parProcIdx_, nsubprobs_proc, MPI::INT,
			subproblemSpecs_, nsubprobsAtRank_, displs, MPI::INT, 0, comm_);

	/** free memory */
	FREE_ARRAY_PTR(displs);

	return STO_RTN_OK;
}

/** solve subproblem */
STO_RTN_CODE DecDdMpi::solveSubproblem()
{
	int doContinue = 1;
	Solutions solutions;
	double stime = MPI::Wtime();

	/** indicators */
	bool doAddFeasCuts  = false;
	bool doAddOptCuts   = false;
	bool doEvalUb       = false;
	bool doCollectSols  = false;
	bool doSolveSubprob = false;
	if (iterCnt_ == 0)
	{
		if (parFeasCuts_ >= 0) doAddFeasCuts = true;
		if (parOptCuts_ >= 0) doAddOptCuts = true;
		if (parEvalUb_ >= 0) doEvalUb = true;
	}
	else
	{
		if (parFeasCuts_ > 0 && iterCnt_ % parFeasCuts_ == 0) doAddFeasCuts = true;
		if (parOptCuts_ > 0 && iterCnt_ % parOptCuts_ == 0) doAddOptCuts = true;
		if (parEvalUb_ > 0 && iterCnt_ % parEvalUb_ == 0) doEvalUb = true;
	}
	if (doAddFeasCuts || doAddOptCuts || doEvalUb)
		doCollectSols = true;
	if (comm_group_ == 0)
		doSolveSubprob = true;

	/**
	 * Solution of subproblem:
	 *  1. solve each subproblem
	 *  2. synchronize solutions
	 *  3. add cover inequalities if the first stage has binary variables only.
	 *  4. generate feasibility cuts; go to 1 if exists.
	 *  5. obtain upper bounds
	 *  6. generate optimality cuts; go to 1 if exists.
	 */

	while (1)
	{
		int doContinueLocal = 1; /**< local flag */

		/** clear solution pool */
		for (unsigned int i = 0; i < solutions.size(); ++i)
			FREE_PTR(solutions[i]);
		solutions.clear();

		/** 1. solve each subproblem */
#ifdef DO_SERIALIZE
		for (int r = 0; r < comm_size_; ++r)
		{
			if (r == comm_rank_)
			{
#endif
				if (doSolveSubprob)
				{
					/** 1. solve each subproblem */
					for (unsigned s = 0; s < subprobs_.size(); ++s)
					{

						if (doContinueLocal == 0) continue;

						if (parFeasCuts_ >= 0 || parOptCuts_ >= 0)
						{
							/** disable cut generator */
							subprobs_[s]->chgCutGenerator(NULL);

							/** update subproblem */
							if (parOptCuts_ >= 0 && primalBound_ < COIN_DBL_MAX)
								subprobs_[s]->updateProblem(NULL, primalBound_);

							/** push cuts */
							subprobs_[s]->pushCuts(cuts_);
							DSPdebugMessage("Rank %d: pushed %d cuts\n", comm_rank_, cuts_->sizeCuts());
						}

						/** update wall clock */
						ticToc();

						/** solution status */
						if (time_remains_ < 0)
						{
							doContinueLocal = 0;
							continue;
						}

						/** set time limit */
						subprobs_[s]->setTimeLimit(CoinMin(time_remains_, parScipTimeLim_));

						/** solve subproblem */
						subprobs_[s]->solve();

						DSPdebugMessage("Rank %d: solved subproblem %d status %d objective %e (%.2f seconds)\n",
								comm_rank_, subprobs_[s]->sind_, subprobs_[s]->si_->getStatus(), subprobs_[s]->si_->getPrimalBound(), MPI::Wtime() - tic_);

						/** update wall clock */
						ticToc();

						/** solution status */
						if (checkStatus(subprobs_[s]) == false || time_remains_ < 0)
						{
							doContinueLocal = 0;
							continue;
						}

						if (doCollectSols)
						{
							/** check duplicate solution */
							CoinPackedVector * x = duplicateSolution(
									model_->getNumSubproblemCouplingCols(s),
									subprobs_[s]->si_->getSolution(),
									ubSolutions_);
							if (x != NULL)
							{
								/** store solution */
								//ubSolutions_.push_back(x);
								solutions.push_back(x);
							}
						}
					}
				}
#ifdef DO_SERIALIZE
			}
			MPI_Barrier(comm_);
		}
#endif

		/** MPI reduce for solution status */
		MPI_Allreduce(&doContinueLocal, &doContinue, 1, MPI::INT, MPI::MIN, comm_);

		/** terminate? */
		if (doContinue == 0)
			break;

		/** 2. synchronize solutions */
		int nsync = 0;
		if (doCollectSols)
			syncSolutions(solutions, nsync);

		/** 3. add cover inequalities if the first stage has "binary variables only" */
		//int nCoverCuts = 0;

		/** 4. generate feasibility cuts; synchronize cuts; go to 1 if violated. */
		if (doAddFeasCuts)
		{
			int ncuts = cuts_->sizeCuts();
			bool violated = false;
			addFeasCuts(solutions);
#ifdef DO_SERIALIZE
			for (int r = 0; r < comm_size_; ++r)
			{
				if (r == comm_rank_)
				{
#endif
					if (doSolveSubprob)
					{
						/** free transform if the current solution is violated. */
						for (unsigned s = 0; s < subprobs_.size(); ++s)
						{
							for (int j = ncuts; j < cuts_->sizeCuts(); ++j)
							{
								double viol = cuts_->rowCutPtr(j)->violated(subprobs_[s]->si_->getSolution());
//								DSPdebugMessage("Rank %d: subproblem %d violation %e\n", comm_rank_, subprobs_[s]->sind_, viol);
								if (viol > 1.0e-6)
								{
//									DSPdebugMessage("Rank %d: free transform %d\n", comm_rank_, subprobs_[s]->sind_);
									violated = true;
									subprobs_[s]->freeTransform();
									break;
								}
							}
						}
					}
#ifdef DO_SERIALIZE
				}
				MPI_Barrier(comm_);
			}
#endif
			if (violated)
				continue;
		}

		/** store feasible solutions */
		for (unsigned int i = 0; i < solutions.size(); ++i)
		{
			ubSolutions_.push_back(solutions[i]);
			solutions[i] = NULL;
		}

		/** 5. obtain upper bounds */
		if (doEvalUb)
		{
			/** update wall clock */
			obtainUpperBounds(nsync);
		}

		/** update wall clock */
		ticToc();

		/** time remains */
		if (time_remains_ < 0)
			doContinueLocal = 0;

		/** MPI reduce for solution status */
		MPI_Allreduce(&doContinueLocal, &doContinue, 1, MPI::INT, MPI::MIN, comm_);

		/** terminate? */
		if (doContinue == 0)
			break;

		/** 6. generate optimality cuts; synchronize cuts; go to 1 if violated. */
		if (doAddOptCuts)
			addOptCuts(nsync);
		break;
	}

	/** clear solution pool */
	for (unsigned int i = 0; i < solutions.size(); ++i)
		FREE_PTR(solutions[i]);
	solutions.clear();

	/** subproblem solution time */
	wtime_subprob_.push_back(MPI::Wtime() - stime);

	/** terminate? */
	if (doContinue == 0)
		return STO_STAT_STOPPED_MPI;

	return STO_RTN_OK;
}

/** synchronize subproblem: send Lagrangian multipliers to subproblems */
STO_RTN_CODE DecDdMpi::syncSubproblem()
{
#define FREE_MEMORY \
	FREE_ARRAY_PTR(sendbuf); \
	FREE_ARRAY_PTR(scounts); \
	FREE_ARRAY_PTR(displs);  \
	FREE_ARRAY_PTR(recvbuf);

	double * sendbuf = NULL;
	int *    scounts = NULL;
	int *    displs  = NULL;
	double * recvbuf = NULL;
	int      rcount = 0;

	if (comm_rank_ == 0)
	{
		int slen = 0; /**< send buffer length */
		for (int s = 0; s < model_->getNumSubproblems(); s++)
			slen += model_->getNumSubproblemCouplingRows(s);

		int curMsglen = 0;

		/** allocate memory */
		sendbuf = new double [slen];
		scounts = new int [comm_size_];
		displs  = new int [comm_size_];

		/** get lambda */
		const double * lambda = master_->getLagrangian();

		/** construct buffer */
		for (int i = 0, scnt = 0; i < comm_size_; ++i)
		{
			scounts[i] = 0;
			for (int j = 0; j < nsubprobsAtRank_[i]; ++j)
			{
				int sind = subproblemSpecs_[scnt];
				//printf("scnt %d sind %d\n", scnt, sind);

				/** split up Lagrangian multipliers for each subproblem */
				const double * msg;
				int msglen = model_->getNumSubproblemCouplingRows(sind); /**< buffer length */
				if (model_->nonanticipativity())
				{
					/** merely to avoid allocating more memory; reconstructing lambda as below should also work */
					msg = lambda + msglen * sind;
				}
				else
				{
					/** construct only the part of lambda that is relevant to the model */
					const int * relevantRows = model_->getSubproblemCouplingRowIndices(sind);
					double * msg_aux = new double [msglen];
					for (int k = 0; k < msglen; k++)
						msg_aux[k] = lambda[relevantRows[k]];
					msg = msg_aux;
				}

#ifdef DSP_DEBUG
				DSPdebugMessage("Sending lambda to rank %d, subproblem %d:  ", i, sind);
				for (int k = 0; k < msglen; k++)
					printf("%.4f ", msg[k]);
				printf("\n");
#endif

				CoinCopyN(msg, msglen, sendbuf + curMsglen);
				curMsglen += msglen;
				scounts[i] += msglen;
				msg = NULL;
				scnt++;
			}
			displs[i] = i == 0 ? 0 : displs[i-1] + scounts[i-1];
		}
		lambda = NULL;
	}
#if USE_ROOT_PROCESSOR == 0
	else
#endif
	{
		rcount = 0;
		for (vector<DecDdSub*>::iterator it = subprobs_.begin(); it != subprobs_.end(); ++it)
		{
			rcount += model_->getNumSubproblemCouplingRows((*it)->sind_);
		}
		recvbuf = new double [rcount];
	}

	/** scatter lambda */
	double stime = MPI::Wtime();
	MPI_Scatterv(sendbuf, scounts, displs, MPI::DOUBLE, recvbuf, rcount, MPI::DOUBLE, 0, comm_);
	DSPdebugMessage("MPI_Scatterv %.2f seconds\n", MPI::Wtime() - stime);

	stime = MPI::Wtime();
#if USE_ROOT_PROCESSOR == 0
	if (comm_rank_ != 0 && comm_group_ == 0)
#else
	if (comm_group_ == 0)
#endif
	{
		int offset = 0;
		int numSubprobs = subprobs_.size();

		for (int i = 0; i < numSubprobs; ++i)
		{
			/** change cut generator */
			if (parFeasCuts_ >= 0 || parOptCuts_ >= 0)
				subprobs_[i]->chgCutGenerator(bdsub_);

			double * lambda = recvbuf + offset;
			offset += model_->getNumSubproblemCouplingRows(subprobs_[i]->sind_);
			subprobs_[i]->updateProblem(lambda, primalBound_);
		}
	}
	DSPdebugMessage("Update subproblems %.2f seconds\n", MPI::Wtime() - stime);

	FREE_MEMORY;

	return STO_RTN_OK;
#undef FREE_MEMORY
}

/** create master problem */
STO_RTN_CODE DecDdMpi::createMaster()
{
	if (comm_rank_ != 0)
		return STO_RTN_OK;

	/** create master problem */
	switch(parMasterAlgo_)
	{
	case Subgradient:
		master_ = new DecDdMasterSubgrad(par_);
		break;
	case DSBM:
		master_ = new DecDdMasterDSB(par_);
		break;
	case IPM:
	case IPM_Feasible:
	default:
		master_ = new DecDdPrimalMaster(par_);
		break;
	}
	master_->createProblem(model_);

	/** allocate memory */
	if (parLogDualVars_)
	{
		multipliers_ = new double [model_->getNumCouplingRows()];
		CoinZeroN(multipliers_, model_->getNumCouplingRows());
	}

	return STO_RTN_OK;
}

/** solve master problem */
STO_RTN_CODE DecDdMpi::solveMaster()
{
	if (comm_rank_ != 0)
		return STO_RTN_OK;

	double stime = MPI_Wtime();

	master_->solve();

	wtime_master_.push_back(MPI_Wtime() - stime);

	if (parLogDualVars_)
	{
		/** calculate change of the Lagrangian multiplier */
		int nmult = model_->getNumCouplingRows();
		double dist = 0.0;
		for (int i = 0; i < nmult; ++i)
		{
			double mult = master_->getLagrangian()[i] - multipliers_[i];
			dist += mult * mult;
		}
		dist = sqrt(dist);
		printf("--> Change of Multiplier %e\n", dist);
		changesOfMultiplier_.push_back(dist);

		/** update multipliers */
		CoinCopyN(master_->getLagrangian(), nmult, multipliers_);
	}

	return STO_RTN_OK;
}

/** synchronize master problem: send subproblem solutions to master and may update dual bound */
STO_RTN_CODE DecDdMpi::syncMaster()
{
#define FREE_MEMORY \
	FREE_ARRAY_PTR(sendbuf); \
	FREE_ARRAY_PTR(recvbuf); \
	FREE_ARRAY_PTR(rcounts); \
	FREE_ARRAY_PTR(displs);

	double * sendbuf = NULL; /**< send buffer */
	double * recvbuf = NULL; /**< receive buffer */
	int *    rcounts = NULL; /**< receive buffer length for each process */
	int *    displs  = NULL; /**< displacement of receive buffer */

	int slenPerProcess = 0;                            /**< send buffer length per process */
	int rbuflen = 0;                                   /**< total length of message buffers from all processes */
	int nsubprobs = model_->getNumSubproblems();
	int * slenPerSubproblem = new int [nsubprobs];
	/** send buffer length per subproblem */
	for (int s = 0; s < nsubprobs; s++)
	{
		slenPerSubproblem[s] = model_->getNumSubproblemCouplingCols(s) + 3;
		rbuflen += slenPerSubproblem[s];
	}

	if (comm_rank_ == 0)
	{
		/** allocate memory */
		recvbuf = new double [rbuflen];
		rcounts = new int [comm_size_];
		displs  = new int [comm_size_];

		/** set displacement */
		displs[0] = 0;
		for (int r = 0, scnt = 0; r < comm_size_; ++r)
		{
			rcounts[r] = 0;
			/** obtain subproblems in each rank */
			for (int j = 0; j < nsubprobsAtRank_[r]; ++j)
			{
				rcounts[r] += slenPerSubproblem[subproblemSpecs_[scnt]];
				scnt++;
			}
			if (r > 0)
				displs[r] = displs[r-1] + rcounts[r-1];
		}
	}
#if USE_ROOT_PROCESSOR == 0
	else
#endif
	if (comm_group_ == 0)
	{
		/** send buffer length per process */
		slenPerProcess = 0;
		for (vector<DecDdSub*>::iterator it = subprobs_.begin(); it != subprobs_.end(); ++it)
			slenPerProcess += slenPerSubproblem[(*it)->sind_];

		/** allocate send buffer memory */
		sendbuf = new double [slenPerProcess];

		/** get message buffer */
		int msgoffset = 0;
		for (vector<DecDdSub*>::iterator it = subprobs_.begin(); it != subprobs_.end(); ++it)
		{
			double * sendbufPtr = sendbuf + msgoffset;
			(*it)->MPImsgbuf(sendbufPtr);
			msgoffset += slenPerSubproblem[(*it)->sind_];
		}
	}

	/** gather subproblem solutions from all processes */
	MPI_Gatherv(sendbuf, slenPerProcess, MPI::DOUBLE,
			recvbuf, rcounts, displs, MPI::DOUBLE, 0, comm_);

	if (comm_rank_ == 0)
	{
		/** allocate memory */
		double * primal_sub   = new double [model_->getNumSubproblems()];
		double * dual_sub   = new double [model_->getNumSubproblems()];
		double ** solution_sub = new double * [model_->getNumSubproblems()];
		for (int s = 0; s < model_->getNumSubproblems(); ++s)
			solution_sub[s] = new double [model_->getNumSubproblemCouplingCols(s)];

		/** parse receive buffer */
		for (int i = 0; i < comm_size_; ++i)
		{
			int offset = 0;
			for (int j = 0; j < nsubprobsAtRank_[i]; ++j)
			{
				/** chop */
				double * msgpiece = recvbuf + displs[i] + offset;

				/** store message */
				int sind = static_cast<int>(msgpiece[0]);
				primal_sub[sind] = msgpiece[1];
				dual_sub[sind] = msgpiece[2];
				CoinCopyN(msgpiece + 3, model_->getNumSubproblemCouplingCols(sind),	solution_sub[sind]);

				offset += slenPerSubproblem[sind];
			}
		}

#ifdef DSP_DEBUG
		DSPdebugMessage("Sending to master:\n");
		for (int i = 0; i < comm_size_; ++i)
		{
			int offset = 0;
			for (int j = 0; j < nsubprobsAtRank_[i]; ++j)
			{
				double * msgpiece = recvbuf + displs[i] + offset;
				int sind = static_cast<int>(msgpiece[0]);
				DSPdebugMessage("s=%d: primobj=%f dualobj=%f, ", sind, primal_sub[sind], dual_sub[sind]);
				for (int k = 0; k < model_->getNumSubproblemCouplingCols(sind); k++)
					printf("%f ", solution_sub[sind][k]);
				printf("\n");

				offset += slenPerSubproblem[sind];
			}
		}
#endif

		/** update master: may update dual bound */
		master_->updateProblem(primalBound_, dualBound_, primal_sub, dual_sub, solution_sub);

		/** store subproblem objective values */
		double primal_subprob = 0.;
		double dual_subprob = 0.;
		for (int s = 0; s < model_->getNumSubproblems(); ++s)
		{
			primal_subprob += primal_sub[s];
			dual_subprob += dual_sub[s];
		}
		primal_subprob_.push_back(primal_subprob);
		dual_subprob_.push_back(dual_subprob);

		/** free memory */
		FREE_ARRAY_PTR(primal_sub);
		FREE_ARRAY_PTR(dual_sub);
		FREE_2D_ARRAY_PTR(model_->getNumSubproblems(), solution_sub);
	}

	FREE_MEMORY;

	return STO_RTN_OK;
#undef FREE_MEMORY
}

/** check whether to terminate or not */
bool DecDdMpi::doTerminate()
{
	bool terminate = false;

	if (comm_rank_ == 0)
	{
		double gapPrimal = 0.0;
		double gapMaster = 0.0;
		gapPrimal = (primalBound_ - dualBound_) / (1.e-10 + fabs(primalBound_));
		if (parMasterAlgo_ != Subgradient)
			gapMaster = (master_->getObjValue() - dualBound_) / (1.e-10 + fabs(master_->getObjValue()));

		/** store primal/dual bounds and master objective value */
		objval_master_.push_back(master_->getObjValue());
		primalBounds_.push_back(primalBound_);
		dualBounds_.push_back(dualBound_);

		/** print iteration information */
		if (parLogLevel_ > 0)
		{
			printf("Iteration %3d: Dual Bound %e", iterCnt_, dualBound_);
			if (parEvalUb_ >= 0)
			{
				printf(", Primal Bound %e ", primalBound_);
				if (gapPrimal > 1.e+4 || gapPrimal < 1.e-4)
					printf("(Optimality Gap %.2e %%)", gapPrimal * 100);
				else
					printf("(Optimality Gap %.2f %%)", gapPrimal * 100);
			}
			printf(", Time Elapsed %.2f sec.", MPI::Wtime() - wtime_start_);
			if (parMasterAlgo_ != Subgradient && parLogLevel_ > 1)
			{
				printf("\n  Master: Obj %e", master_->getObjValue());
				if (parLogLevel_ > 2)
				{
					if (gapMaster > 1.e+4 || gapMaster < 1.e-4)
						printf(" (Gap %.2e %%)", gapMaster * 100);
					else
						printf(" (Gap %.2f %%)", gapMaster * 100);
				}
				printf(", Solution Time %.2f sec.", wtime_master_[wtime_master_.size() - 1]);
			}
			printf("\n");
		}

		/** If master is not optimal, terminate. */
		status_ = master_->getStatus();

		if (parIterLim_ <= iterCnt_)
		{
			status_ = STO_STAT_LIM_ITERorTIME;
			terminate = true;
		}
		else if (status_ != STO_STAT_OPTIMAL)
			terminate = true;
		else if (gapPrimal < parStopTol_)
			terminate = true;
		else
			terminate = master_->terminate();
	}

	/** broadcast decision */
	MPI_Bcast(&terminate, 1, MPI::BOOL, 0, comm_);

	return terminate;
}

/** collect results */
STO_RTN_CODE DecDdMpi::collectResults()
{
	/** number of iterations */
	numIterations_ = iterCnt_;

	/** solution time */
	solutionTime_ = MPI_Wtime() - wtime_start_;

	/** cpu time */
	ctime_start_ = CoinCpuTime() - ctime_start_;

	/** gather cpu time */
	double total_ctime = 0.;
	MPI_Allreduce(&ctime_start_, &total_ctime, 1, MPI_DOUBLE, MPI_SUM, comm_);
	ctime_start_ = total_ctime;

#ifdef TSSDD_WRITE_FILE
	char logfile[128];

	/** write results in file */
	sprintf(logfile, "%sDecDdMpi%d.log", logfile_prefix_.c_str(), comm_rank_);

	ofstream f(logfile, ofstream::out);
#endif

	if (comm_rank_ == 0)
	{
#ifdef TSSDD_WRITE_FILE
		/** write results in file */
		char logfileMaster[128];
		sprintf(logfileMaster, "%sDecDdMpiMaster.log", logfile_prefix_.c_str());

		int niters = master_->wtime_solve_.size();

		ofstream f2(logfileMaster, ofstream::out);
		f2 << " \tIter\tTimeElapsed\tMasterObjval\tDualObjval\tPrimalBound\tDualBound\tSolutionTime" << endl;
		for (int i = 0; i < niters; ++i)
		{
			f2 << i << "\t"
					<< master_->itncnt_[i] << "\t"
					<< wtime_elapsed_[i] << "\t"
					<< master_->objvals_[i] << "\t"
					<< master_->dual_objvals_[i] << "\t"
					<< master_->primal_bounds_[i] << "\t"
					<< master_->dual_bounds_[i] << "\t"
					<< master_->wtime_solve_[i] << endl;
		}
		f2.close();
#endif
	}
#ifdef TSSDD_WRITE_FILE
#if USE_ROOT_PROCESSOR == 0
	else
#endif
	{
		int nsubprobs = subprobs_.size();
		if (nsubprobs > 0)
		{
#if 0
			int niters = subprobs_[0]->wtime_solution_.size();
			f << "Iter";
			for (int j = 0; j < nsubprobs; ++j)
				f << "\tScenario " << subprobs_[j]->sind_;
			f << endl;
			for (int i = 0; i < niters; ++i)
			{
				f << i;
				for (int j = 0; j < nsubprobs; ++j)
					f << "\t" << subprobs_[j]->wtime_solution_[i];
				f << endl;
			}
#else
			int niters = wtime_subprob_.size();
			f << "Iter\tTime" << endl;
			for (int i = 0; i < niters; ++i)
				f << i << "\t" << wtime_subprob_[i] << endl;
#endif
		}
	}

	f.close();
#endif

	return STO_RTN_OK;
}

/** get upper bound */
double DecDdMpi::getUpperBound(
		DecDdSub * ddsub, /**< subproblem to evaluate */
		bool & feasible   /**< indicating feasibility */)
{
	assert(bdsub_);

	/** get solution */
	const double * solution = ddsub->si_->getSolution();

	return getUpperBound(solution, feasible);
}

/** get upper bound */
double DecDdMpi::getUpperBound(
		const double * solution, /**< solution to evaluate */
		bool & feasible          /**< indicating feasibility */)
{
	assert(solution);

	StoModel * stoModel = dynamic_cast<StoModel *>(model_);
	assert(stoModel != NULL); /** this function should only be called for stochastic models */

	/** get objective value for the first stage */
	double objval = 0.0;
	double objval0 = 0.0;
	for (int j = 0; j < stoModel->getNumCols(0); ++j)
		objval0 += stoModel->getObjCore(0)[j] * solution[j];

	/** get recourse value */
	feasible = true;
	int nsubprobs = subprobs_.size();
	double objval_reco = COIN_DBL_MAX;
	for (int s = 0; s < nsubprobs; ++s)
	{
		/** solve recourse problem */
		bdsub_->solveSingleRecourse(subprobs_[s]->sind_, solution, objval_reco, NULL);

		/** infeasible ?*/
		if (bdsub_->status_[subprobs_[s]->sind_] == STO_STAT_PRIM_INFEASIBLE)
		{
			feasible = false;
			break;
		}
		/** objective */
		objval += (objval0 + objval_reco) * stoModel->getProbability()[subprobs_[s]->sind_];
		DSPdebugMessage("Rank %d scenario %d objective(first) %e objective(reco) %e objective(wsum) %e\n",
				comm_rank_, subprobs_[s]->sind_, objval0, objval_reco, (objval0 + objval_reco) * stoModel->getProbability()[subprobs_[s]->sind_]);
	}

	/** infeasible? */
	if (!feasible)
		return COIN_DBL_MAX;

	return objval;
}

/** check solution status and determine whether to continue or stop. */
bool DecDdMpi::checkStatus(DecDdSub * ddsub)
{
	bool doContinue = false;
	switch (ddsub->si_->getStatus())
	{
	case STO_STAT_OPTIMAL:
	case STO_STAT_LIM_ITERorTIME:
	case STO_STAT_STOPPED_GAP:
	case STO_STAT_STOPPED_NODE:
	case STO_STAT_STOPPED_TIME:
		doContinue = true;
		break;
	default:
		printf("Warning: solution status is %d\n", ddsub->si_->getStatus());
		break;
	}
	return doContinue;
}

/** check whether solution is duplicate or not; return NULL if duplicate */
CoinPackedVector * DecDdMpi::duplicateSolution(
		int size,           /**< size of array */
		const double * x,   /**< current solution */
		Solutions solutions /**< solution pool to check duplication */)
{
	assert(x);
	bool dup = false;

	CoinPackedVector * xvec = new CoinPackedVector;
	for (int i = 0; i < size; ++i)
	{
		if (fabs(x[i]) > 1.0e-8)
			xvec->insert(i, x[i]);
	}
	dup = duplicateVector(xvec, solutions);

	/** free if duplicate */
	if (dup) FREE_PTR(xvec);

	return xvec;
}

/** synchronize solution pool */
void DecDdMpi::syncSolutions(
		Solutions & solutions,
		int & numSyncedSolutions)
{
#ifdef DSP_DEBUG
	int numBeforeSync = solutions.size();
#endif
	Solutions solutionsToSync;
	Solutions solutionsGathered;

	/** vector of solutions to sync */
	for (unsigned i = numSyncedSolutions; i < solutions.size(); ++i)
		solutionsToSync.push_back(solutions[i]);

	/** erase vector elements moved to solutionsToSync*/
	solutions.erase(solutions.begin() + numSyncedSolutions, solutions.end());

	/** After this call, every processor should have the same solutionsGathered. */
	MPIscatterCoinPackedVectors(comm_, solutionsToSync, solutionsGathered);

	/** move solutions */
	for (unsigned i = 0; i < solutionsGathered.size(); ++i)
		solutions.push_back(solutionsGathered[i]);

#if 0
	/** print solution characteristics */
	for (int r = 0; r < comm_size_; ++r)
	{
		if (r == comm_rank_)
		{
			for (unsigned i = 0; i < solutionsGathered.size(); ++i)
			{
				printf("Rank %d: solution[%d] nzcnt %d minindex %d maxindex %d sum %e 1-norm %e 2-norm %e inf-norm %e\n",
						comm_rank_, i, solutionsGathered[i]->getNumElements(),
						solutionsGathered[i]->getMinIndex(), solutionsGathered[i]->getMaxIndex(),
						solutionsGathered[i]->sum(), solutionsGathered[i]->oneNorm(),
						solutionsGathered[i]->twoNorm(), solutionsGathered[i]->infNorm());
			}
		}
		MPI_Barrier(comm_);
	}
#endif

	/** update number of solutions synchronized */
	numSyncedSolutions = solutions.size();

	/** delete local solutions */
	for (unsigned i = 0; i < solutionsToSync.size(); ++i)
		FREE_PTR(solutionsToSync[i]);

#if 0
	int nSolutions = solutionsGathered.size();
	for (int i = 0; i < nSolutions; ++i)
	{
		if (duplicateVector(solutionsGathered[i], solutions) == false)
			solutions.push_back(solutionsGathered[i]);
		else
			FREE_PTR(solutionsGathered[i]);
	}
	numSyncedSolutions = solutions.size();
#endif

#ifdef DSP_DEBUG
	int numAfterSync = solutions.size();
	DSPdebugMessage("Rank %d: received %d new solutions from the other processes. total solutions %d.\n",
			comm_rank_, numAfterSync - numBeforeSync, numAfterSync);
#endif
}

/** obtain upper bounds for solution pools */
void DecDdMpi::obtainUpperBounds(
		int num /**< number of solutions to evaluate */)
{
	StoModel * stoModel = dynamic_cast<StoModel *>(model_);
	assert(stoModel != NULL); /** this function should only be called for stochastic models */

	/** 1. Evaluate the solutions; */
	int timeExceeded = 0;
	int timeExceededLocal = 0;
	bool feasible = true;
	double * ub = new double [num];
	CoinZeroN(ub, num);

#ifdef DO_SERIALIZE
	for (int r = 0; r < comm_size_; r++)
	{
		if (r == comm_rank_)
		{
#endif
			for (int i = comm_group_; i < num; i += num_comm_groups_)
			{
				double * sol = ubSolutions_[numSyncedUbSolutions_ + i]->denseVector(stoModel->getNumCols(0));
				ub[i] = getUpperBound(sol, feasible);
				DSPdebugMessage("Rank %d: ub[%d] %e %s\n", comm_rank_, i, ub[i], feasible ? "feas" : "infeas");
				FREE_ARRAY_PTR(sol);

				ticToc();
				if (time_remains_ < 0.0)
				{
					timeExceededLocal = 1;
					break;
				}
			}
#ifdef DO_SERIALIZE
		}
		MPI_Barrier(comm_);
	}
#endif
	
	MPI_Allreduce(&timeExceededLocal, &timeExceeded, 1, MPI::INT, MPI::MAX, comm_);
	if (timeExceeded == 1)
	{
		FREE_ARRAY_PTR(ub);
		return;
	}

	/** 2. Collect the objective values; */
	collectObjValues(num, ub);

	/** 3. Take the minimum */
	bool ubUpdated = false;
	for (int i = 0; i < num; ++i)
	{
		DSPdebugMessage("Rank %d: ub[%d] %e\n", comm_rank_, i, ub[i]);
		if (primalBound_ > ub[i])
		{
			primalBound_ = ub[i];
			ubUpdated = true;
		}
	}
	FREE_ARRAY_PTR(ub);
	numSyncedUbSolutions_ += num;

	/** 4. print */
	if (parLogLevel_ > 2)
	{
		if (ubUpdated)
			printf(" -> BEST Primal Bound %e\n", primalBound_);
	}
	else if (parLogLevel_ == 1)
	{
		if (ubUpdated) printf(" * ");
		else printf("   ");
	}
}

/** collect objective values;
 * This takes the summation of all vals from different processors. */
void DecDdMpi::collectObjValues(
		int num,      /**< size of vals */
		double * vals /**< objective values */)
{
	double * sumvals = new double [num];

	/** MPI_Allreduce */
	MPI_Allreduce(vals, sumvals, num, MPI::DOUBLE, MPI::SUM, comm_);

	/** copy values */
	CoinCopyN(sumvals, num, vals);

	/** free memory */
	FREE_ARRAY_PTR(sumvals);
}

/** add feasibility cuts */
int DecDdMpi::addFeasCuts(
		Solutions solutions /**< solutions to check */)
{
	StoModel * stoModel = dynamic_cast<StoModel *>(model_);
	assert(stoModel != NULL); /** this function should only be called for stochastic models */

	int ncols_first = stoModel->getNumCols(0); /**< number of first-stage variables */
	int nsols = solutions.size();            /**< number of solutions stored */
	double ** sol = new double * [nsols];    /**< array of dense vector of solutions */
	OsiCuts * cuts = new OsiCuts;            /**< cuts generated */

	/** loop over solutions */
	for (int i = 0; i < nsols; ++i)
		sol[i] = NULL;
#ifdef DO_SERIALIZE
	for (int r = 0; r < comm_size_; ++r)
	{
		if (r == comm_rank_)
		{
			DSPdebugMessage("Rank %d generating feasibility cuts\n", comm_rank_);
#endif
			for (int i = comm_group_; i < nsols; i += num_comm_groups_)
			{
				sol[i] = solutions[i]->denseVector(ncols_first);
				bdsub_->generateCuts(ncols_first, sol[i], cuts, TssBdSub::TssDd, TssBdSub::FeasCut);
			}
#ifdef DO_SERIALIZE
		}
		MPI_Barrier(comm_);
	}
#endif


#ifdef DSP_DEBUG
	/** print cuts */
	DSPdebugMessage("Rank %d: found %d feasibility cuts\n", comm_rank_, cuts->sizeCuts());
	for (int i = 0; i < cuts->sizeCuts(); ++i)
	{
		OsiRowCut * rc = cuts->rowCutPtr(i);
		CoinPackedVector rcrow = rc->row();
		rc->print();
		DSPdebugMessage("Rank %d: cut[%d] rhs %e minind %d maxind %d nzcnt %d sum %e 1-norm %e 2-norm %e inf-norm %e\n",
				comm_rank_, i, rc->lb(), rcrow.getMinIndex(), rcrow.getMaxIndex(), rcrow.getNumElements(), rcrow.sum(), rcrow.oneNorm(), rcrow.twoNorm(), rcrow.infNorm());
	}
#endif

	/** synchronize cuts */
	syncCuts(0, cuts);

	/** add generated cuts to pool */
	for (int i = 0; i < cuts->sizeCuts(); ++i)
		cuts_->insertIfNotDuplicate(cuts->rowCut(i));
	numSyncedCuts_ = cuts_->sizeCuts();

	DSPdebugMessage("Rank %d: feasibility cut synchronized %d, cut pool %d\n", comm_rank_, cuts->sizeCuts(), numSyncedCuts_);

	/** free memory */
	FREE_PTR(cuts);
	FREE_2D_ARRAY_PTR(nsols, sol);

	return 1;
}

/** add optimality cuts */
int DecDdMpi::addOptCuts(
		int num /**< number of solutions to check */)
{
	StoModel * stoModel = dynamic_cast<StoModel *>(model_);
	assert(stoModel != NULL); /** this function should only be called for stochastic models */

	int violated = 0;
	int ncols_first = stoModel->getNumCols(0);    /**< number of first-stage variables */
	int nSolutions = ubSolutions_.size();       /**< number of solutions stored */
	vector<int> isol;                           /**< solution indices considered in process */
	int * numsols   = new int [comm_size_];     /**< number of solutions considered in each process */
	int * displs    = new int [comm_size_];
	int * cutgroup  = NULL;
	double ** sol   = new double * [num];       /**< array of dense vector of solutions */
	OsiCuts * cuts  = new OsiCuts;              /**< cuts generated */
	double * aggrow = new double [ncols_first]; /**< aggregated optimality cut row */
	double aggrhs = 0.0;                        /**< aggregated optimality cut rhs */

	/** loop over solutions */
	for (int i = 0; i < num; ++i)
		sol[i] = NULL;
#ifdef DO_SERIALIZE
	for (int r = 0; r < comm_size_; ++r)
	{
		if (r == comm_rank_)
		{
			DSPdebugMessage("Rank %d generating optimality cuts\n", comm_rank_);
#endif
			for (int i = comm_group_; i < num; i += num_comm_groups_)
			{
				isol.push_back(i);
				sol[i] = ubSolutions_[nSolutions - num + i]->denseVector(ncols_first);
				bdsub_->generateCuts(ncols_first + 1, sol[i], cuts, TssBdSub::TssDd, TssBdSub::OptCut);
			}
#ifdef DO_SERIALIZE
		}
		MPI_Barrier(comm_);
	}
#endif

#ifdef DSP_DEBUG
	/** print cuts */
	for (int r = 0; r < comm_size_; ++r)
	{
		if (r == comm_rank_)
		{
			DSPdebugMessage("Rank %d: found %d optimality cuts\n", comm_rank_, cuts->sizeCuts());
			for (int i = 0; i < cuts->sizeCuts(); ++i)
			{
				OsiRowCut * rc = cuts->rowCutPtr(i);
				CoinPackedVector rcrow = rc->row();
				rc->print();
				DSPdebugMessage("Rank %d: cut[%d] rhs %e minind %d maxind %d nzcnt %d sum %e 1-norm %e 2-norm %e inf-norm %e\n",
						comm_rank_, i, rc->lb(), rcrow.getMinIndex(), rcrow.getMaxIndex(), rcrow.getNumElements(), rcrow.sum(), rcrow.oneNorm(), rcrow.twoNorm(), rcrow.infNorm());
			}
		}
		MPI_Barrier(comm_);
	}
#endif

	/** synchronize cuts */
	syncCuts(0, cuts);

	/** synchronize solution indices */
	int isol_size[1] = {(int) isol.size()};
	MPI_Allgather(isol_size, 1, MPI::INT, numsols, 1, MPI::INT, comm_);

	for (int i = 0; i < comm_size_; ++i)
		displs[i] = i == 0 ? 0 : displs[i-1] + numsols[i-1];
	cutgroup = new int [cuts->sizeCuts()];
	MPI_Allgatherv(&isol[0], (int) isol.size(), MPI::INT, cutgroup, numsols, displs, MPI::INT, comm_);

	/** construct cuts */
	for (int k = 0; k < num; ++k)
	{
		/** initializing */
		for (int j = 0; j < ncols_first; ++j)
			aggrow[j] = -(stoModel->getObjCore(0)[j]);
		aggrhs = 0.0;

		/** aggregating */
		for (int i = 0; i < cuts->sizeCuts(); ++i)
		{
			if (cutgroup[i] != k) continue;

			OsiRowCut * rc = cuts->rowCutPtr(i);
			if (rc == NULL) continue;

			/** retrieve cut data */
			CoinPackedVector row = rc->row();

			/** aggregate cut rows */
			for (int j = 0; j < row.getNumElements(); ++j)
			{
				if (row.getIndices()[j] < ncols_first)
					aggrow[row.getIndices()[j]] += row.getElements()[j];
//				else
//					printf("Warning: cut row index exceeds %d\n", row.getIndices()[j]);
			}

			/** aggregate cut rhs */
			aggrhs += rc->lb();
		}

		/** create new cut */
		CoinPackedVector cutrow;
		OsiRowCut rowcut;

		for (int j = 0; j < ncols_first; ++j)
		{
			if (fabs(aggrow[j]) > 1E-10)
				cutrow.insert(j, aggrow[j]);
		}
		cutrow.insert(ncols_first + stoModel->getNumCols(1), 1.0);

		rowcut.setRow(cutrow);
		rowcut.setUb(COIN_DBL_MAX);
		rowcut.setLb(aggrhs);
		//rowcut.print();

		/** store cut */
		cuts_->insertIfNotDuplicate(rowcut);
	}
	numSyncedCuts_ = cuts_->sizeCuts();

	DSPdebugMessage("Rank %d: optimality cut synchronized %d, cut pool %d\n", comm_rank_, cuts->sizeCuts(), numSyncedCuts_);

	/** check the violation */
#if 0
	for (int i = 0; i < num; ++i)
	{
		for (int j = nCuts; j < cuts_->sizeCuts(); ++j)
		{
			double eff = cuts_->rowCutPtr(j)->violated(sol[i]) / cuts_->rowCutPtr(j)->row().twoNorm();
			if (eff > 1.0e-6)
			{
				violated = 1;
				break;
			}
		}
		if (violated) break;
	}
#endif

	/** free memory */
	FREE_PTR(cuts);
	FREE_PTR(numsols);
	FREE_PTR(displs);
	FREE_PTR(cutgroup);
	FREE_2D_ARRAY_PTR(num, sol);
	FREE_ARRAY_PTR(aggrow);

	return violated;
}

/** synchronize cut pool */
void DecDdMpi::syncCuts(
		int start,     /**< start index to synchronize */
		OsiCuts * cuts /**< cuts to synchronize */)
{
	OsiCuts cutsToSync;
	OsiCuts cutsGathered;

	/** cuts to sync */
	//DSPdebugMessage("Rank %d: number of cuts %d before synchronization\n", comm_rank_, (int) cuts_->sizeCuts());
	for (int i = start; i < cuts->sizeCuts(); ++i)
	{
		OsiRowCut * rc = cuts->rowCutPtr(i);
		cutsToSync.insert(rc);
	}

	/** clean-up cut pool */
	OsiCuts tmpcuts;
	for (int i = 0; i < start; ++i)
	{
		OsiRowCut * rc = cuts->rowCutPtr(i);
		tmpcuts.insert(rc);
	}
	cuts->dumpCuts();
	for (int i = 0; i < start; ++i)
	{
		OsiRowCut * rc = tmpcuts.rowCutPtr(i);
		cuts->insert(rc);
	}
	tmpcuts.dumpCuts();

	/** synchronize cuts */
	MPIscatterOsiCuts(comm_, cutsToSync, cutsGathered);

	/** move cuts gathered */
	for (int i = 0; i < cutsGathered.sizeCuts(); ++i)
	{
		OsiRowCut * rc = cutsGathered.rowCutPtr(i);
		cuts->insert(rc);
	}

	/** dump cuts so the original cut pointers were not released. */
	cutsGathered.dumpCuts();
}

/** update wall clock time remains */
void DecDdMpi::ticToc()
{
	time_remains_ -= MPI::Wtime() - tic_;
	tic_ = MPI::Wtime();
}

