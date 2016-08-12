/*
 * TssDe.cpp
 *
 *  Created on: Sep 24, 2014
 *      Author: kibaekkim
 */

//#define DSP_DEBUG
//#define DSP_TIMING

/** DSP */
#include "Utility/StoMessage.h"
#include "Solver/TssDe.h"
#include "Solver/SolverInterfaceClp.h"
#include "Solver/SolverInterfaceScip.h"

TssDe::TssDe():
TssSolver(),
si_(NULL)
{
	/** nothing to do */
}

TssDe::~TssDe()
{
	FREE_PTR(si_)
}

/** solve */
STO_RTN_CODE TssDe::solve()
{
#define FREE_MEMORY       \
	FREE_PTR(mat)         \
	FREE_ARRAY_PTR(clbd)  \
	FREE_ARRAY_PTR(cubd)  \
	FREE_ARRAY_PTR(ctype) \
	FREE_ARRAY_PTR(obj)   \
	FREE_ARRAY_PTR(rlbd)  \
	FREE_ARRAY_PTR(rubd)

	assert(model_);

	/** model info */
	CoinPackedMatrix * mat = NULL;
	double * clbd   = NULL;
	double * cubd   = NULL;
	double * obj    = NULL;
	char *   ctype  = NULL;
	double * rlbd   = NULL;
	double * rubd   = NULL;

	BGN_TRY_CATCH

	double stime;

#ifdef DSP_TIMING
	stime = CoinGetTimeOfDay();
#endif
	/** get De model */
	STO_RTN_CHECK_THROW(
			model_->copyDeterministicEquivalent(mat, clbd, cubd, ctype, obj, rlbd, rubd),
			"copyDeterministicEquivalent", "TssModel");
#ifdef DSP_TIMING
	printf("Constructed extensive form data, elapsed time %.2f seconds\n", CoinGetTimeOfDay() - stime);
	stime = CoinGetTimeOfDay();
#endif

	int nIntegers = model_->getNumCoreIntegers();

	/** relax integrality? */
	if (par_->relaxIntegrality_[0])
	{
		for (int j = 0; j < model_->getNumCols(0); ++j)
		{
			if (ctype[j] != 'C')
				nIntegers--;
			ctype[j] = 'C';
		}
	}
	if (par_->relaxIntegrality_[1])
	{
		for (int j = 0; j < model_->getNumCols(1); ++j)
		{
			if (ctype[model_->getNumCols(0) + j] != 'C')
				nIntegers--;
		}
		CoinFillN(ctype + model_->getNumCols(0), model_->getNumScenarios() * model_->getNumCols(1), 'C');
	}

	if (nIntegers > 0)
	{
		si_ = new SolverInterfaceScip(par_);
		/** print level */
		si_->setPrintLevel(CoinMin(par_->logLevel_ + 2, 5));
	}
	else
	{
		si_ = new SolverInterfaceClp(par_);
		/** print level */
		si_->setPrintLevel(par_->logLevel_);
	}

	/** load problem */
	si_->loadProblem(mat, clbd, cubd, obj, ctype, rlbd, rubd, "TssDe");
#ifdef DSP_TIMING
	printf("Loaded problem to solver, elapsed time %.2f seconds\n", CoinGetTimeOfDay() - stime);
#endif

	/** time limit */
	si_->setTimeLimit(par_->wtimeLimit_);

	/** set iteration limit */
	si_->setIterLimit(par_->iterLimit_);

	/** set node limit */
	si_->setNodeLimit(par_->nodeLimit_);

	/** solution time */
	stime = clockType_ ? CoinGetTimeOfDay() : CoinCpuTime();

	/** solve */
	si_->solve();

	/** solution time */
	solutionTime_ = (clockType_ ? CoinGetTimeOfDay() : CoinCpuTime()) - stime;

	/** solution status */
	status_ = si_->getStatus();
printf("status %d\n", status_);
	/** get solutions */
	if (status_ == STO_STAT_OPTIMAL ||
		status_ == STO_STAT_STOPPED_TIME ||
		status_ == STO_STAT_STOPPED_NODE ||
		status_ == STO_STAT_STOPPED_GAP)
	{
		/** objective bounds */
		primalBound_ = si_->getPrimalBound();
		dualBound_ = si_->getDualBound();

		/** solution */
		assert(solution_);
		if (si_->getSolution())
			CoinCopyN(si_->getSolution(), si_->getNumCols(), solution_);

		/** statistics */
		numIterations_ = si_->getIterationCount();
		numNodes_ = si_->getNumNodes();
	}

	/** save memory */
	FREE_MEMORY

	END_TRY_CATCH_RTN(FREE_MEMORY,STO_RTN_ERR)

	return STO_RTN_OK;

#undef FREE_MEMORY
}

void TssDe::writeMps(const char* filename)
{
#define FREE_MEMORY       \
	FREE_PTR(mat)         \
	FREE_ARRAY_PTR(clbd)  \
	FREE_ARRAY_PTR(cubd)  \
	FREE_ARRAY_PTR(ctype) \
	FREE_ARRAY_PTR(obj)   \
	FREE_ARRAY_PTR(rlbd)  \
	FREE_ARRAY_PTR(rubd)

	assert(model_);

	/** model info */
	CoinPackedMatrix * mat = NULL;
	double * clbd   = NULL;
	double * cubd   = NULL;
	double * obj    = NULL;
	char *   ctype  = NULL;
	double * rlbd   = NULL;
	double * rubd   = NULL;

	BGN_TRY_CATCH

	double stime;

	/** get De model */
	STO_RTN_CHECK_THROW(
			model_->copyDeterministicEquivalent(mat, clbd, cubd, ctype, obj, rlbd, rubd),
			"copyDeterministicEquivalent", "TssModel");

	int nIntegers = model_->getNumCoreIntegers();

	if (nIntegers > 0)
	{
		si_ = new SolverInterfaceScip(par_);
		/** print level */
		si_->setPrintLevel(CoinMin(par_->logLevel_ + 2, 5));
	}
	else
	{
		si_ = new SolverInterfaceClp(par_);
		/** print level */
		si_->setPrintLevel(par_->logLevel_);
	}

	/** load problem */
	si_->loadProblem(mat, clbd, cubd, obj, ctype, rlbd, rubd, "TssDe");

	/** write Mps file  */
	si_->writeMps(filename);

	/** save memory */
	FREE_MEMORY

	END_TRY_CATCH_RTN(FREE_MEMORY,;)

#undef FREE_MEMORY
}
