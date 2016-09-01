/*
 * TssDdPrimalMaster.cpp
 *
 *  Created on: Dec 11, 2014
 *      Author: kibaekkim
 */

/** DSP */
#include "Utility/StoMessage.h"
#include "Utility/StoUtility.h"
#include "Solver/TssDdPrimalMaster.h"
#include "Solver/SolverInterfaceSpx.h"
#include "Solver/SolverInterfaceClp.h"
#include "Solver/OoqpEps.h"

#include "CoinWarmStartBasis.hpp"

#define ENABLE_UPPER_BOUNDING 0

TssDdPrimalMaster::~TssDdPrimalMaster()
{
	FREE_PTR(si_);
	FREE_PTR(cuts_);
	FREE_ARRAY_PTR(prox_);
	FREE_ARRAY_PTR(solution_);
	/** free cuts */
	for (int i = 0; i < ubCuts_->sizeCuts(); ++i)
		delete ubCuts_->rowCutPtr(i);
	ubCuts_->dumpCuts();
	FREE_PTR(ubCuts_);
}

/** create problem */
STO_RTN_CODE TssDdPrimalMaster::createProblem(const TssModel * model)
{
#define FREE_MEMORY       \
	FREE_ARRAY_PTR(ctype); \
	FREE_ARRAY_PTR(clbd); \
	FREE_ARRAY_PTR(cubd); \
	FREE_ARRAY_PTR(obj);  \
	FREE_ARRAY_PTR(rlbd); \
	FREE_ARRAY_PTR(rubd); \
	FREE_ARRAY_PTR(bgn);  \
	FREE_ARRAY_PTR(len);  \
	FREE_ARRAY_PTR(ind);  \
	FREE_ARRAY_PTR(elem); \
	FREE_PTR(mat);

	assert(model);

	int i, nzcnt, pos;
	char * ctype = NULL;
	double * clbd = NULL;
	double * cubd = NULL;
	double * obj  = NULL;
	double * rlbd = NULL;
	double * rubd = NULL;
	CoinBigIndex * bgn = NULL;
	int * len = NULL;
	int * ind = NULL;
	double * elem = NULL;
	CoinPackedMatrix * mat = NULL;

	BGN_TRY_CATCH

	/** retrieve model informatio */
	nscen_ = model->getNumScenarios();
	ncols_first_ = model->getNumCols(0);
	nCutsPerIter_ = CoinMin(nCutsPerIter_, nscen_);

	/** LP dimension */
	nrows_ = ncols_first_;
	ncols_ = nCutsPerIter_ + nscen_ * ncols_first_;
#if ENABLE_UPPER_BOUNDING == 1
	if (par_->TssDdFeasRecoveryInt_ >= 0)
		ncols_++;
#endif
	nzcnt  = nscen_ * ncols_first_;

	/** allocate memory */
	ctype = new char [ncols_];
	clbd = new double[ncols_];
	cubd = new double[ncols_];
	obj  = new double[ncols_];
	rlbd = new double[nrows_];
	rubd = new double[nrows_];
	bgn  = new CoinBigIndex[nrows_ + 1];
	len  = new int[nrows_];
	ind  = new int[nzcnt];
	elem = new double[nzcnt];

	/** all continuous variables */
	CoinFillN(ctype, ncols_, 'C');

	/** c */
	CoinFillN(obj, nCutsPerIter_, 1.0);
	CoinZeroN(obj + nCutsPerIter_, nscen_ * ncols_first_);
#if ENABLE_UPPER_BOUNDING == 1
	if (par_->TssDdFeasRecoveryInt_ >= 0)
		obj[ncols_ - 1] = 0.0;
#endif

	/** proximal point for trust region */
	if (par_->TssDdEnableTrustRegion_)
	{
		assert(rho_ < COIN_DBL_MAX);
		prox_ = new double [nscen_ * ncols_first_];
		CoinZeroN(prox_, nscen_ * ncols_first_);
	}
	else
		rho_ = COIN_DBL_MAX;

	/** bounds */
	CoinFillN(clbd, nCutsPerIter_, -COIN_DBL_MAX);
	CoinFillN(cubd, nCutsPerIter_, +COIN_DBL_MAX);

	/** trust region */
	CoinFillN(clbd + nCutsPerIter_, nscen_ * ncols_first_, -rho_);
	CoinFillN(cubd + nCutsPerIter_, nscen_ * ncols_first_, rho_);

#if ENABLE_UPPER_BOUNDING == 1
	if (par_->TssDdFeasRecoveryInt_ >= 0)
	{
		clbd[ncols_ - 1] = -COIN_DBL_MAX;
		cubd[ncols_ - 1] = COIN_DBL_MAX;
	}
#endif

	/** row bounds */
	CoinZeroN(rlbd, nrows_);
	CoinZeroN(rubd, nrows_);

	/** for constraints */
	pos = 0;
	for (i = 0; i < ncols_first_; ++i)
	{
		bgn[i] = pos;
		for (int j = 0; j < nscen_; ++j)
		{
			ind[pos] = nCutsPerIter_ + j * ncols_first_ + i;
			elem[pos] = 1.0;
			pos++;
		}
		len[i] = pos - bgn[i];
	}
	bgn[nrows_] = pos;
	assert(pos == nzcnt);

	/** constraint matrix */
	mat = new CoinPackedMatrix(false, ncols_, nrows_, nzcnt, elem, ind, bgn, len);
	DSPdebug(mat->verifyMtx(4));

	/** create solver interface */
	switch (par_->TssDdMasterSolver_)
	{
	case Simplex:
		si_ = new SolverInterfaceClp(par_);
		break;
	case IPM:
		si_ = new SolverInterfaceOoqp(par_);
		break;
	case IPM_Feasible:
		si_ = new OoqpEps(par_);
		break;
	default:
		si_ = new OoqpEps(par_);
		break;
	}

	/** maximization */
	si_->setObjSense(-1);

	/** copy problem data */
	si_->loadProblem(mat, clbd, cubd, obj, ctype, rlbd, rubd);

	/** allocate memory for solution */
	solution_ = new double [ncols_];
	CoinZeroN(solution_, ncols_);

	/** initialize cut pool */
	cuts_ = new OsiCuts;
	ubCuts_ = new OsiCuts;

	/** set print level */
	setPrintLevel(CoinMax(0, par_->logLevel_ - 2));

	END_TRY_CATCH_RTN(FREE_MEMORY,STO_RTN_ERR)

	FREE_MEMORY;

	return STO_RTN_OK;
#undef FREE_MEMORY
}

/** update problem: may update dual bound */
STO_RTN_CODE TssDdPrimalMaster::updateProblem(
		double primal_bound, /**< primal bound of the original problem */
		double & dual_bound, /**< dual bound of the original problem */
		double * objvals,    /**< objective values of subproblems */
		double ** solution   /**< subproblem solutions */)
{
//#define DEBUG_TRUST_REGION

	int nCutsAdded = 0;
	double primalBound = si_->getPrimalBound();

	/** aging cuts */
	agingCuts();

	double newobj = 0.0;
	for (int s = 0; s < nscen_; ++s)
	{
#if 0
		printf("Scenario %d objective %e\n", s, objvals[s]);
		for (int j = 0, jj = 0; j < ncols_first_; ++j)
		{
			if (jj > 0 && jj % 5 == 0) printf("\n");
			if (fabs(solution[s][j]) > 1.e-10)
			{
				printf("%4d [%e] ", j, solution[s][j]);
				jj++;
			}
		}
		printf("\n");
#endif
		newobj += objvals[s];
	}

	/** set dual objective value */
	dual_objvals_.push_back(newobj);

#if 0
	/** print proximal point */
	printf("proximal point:\n");
	for (int j = 0, jj = 0; j < nscen_ * ncols_first_; ++j)
	{
		if (jj > 0 && jj % 5 == 0) printf("\n");
		if (fabs(prox_[j]) > 1.e-10)
		{
			printf("%4d [%e] ", j, prox_[j]);
			jj++;
		}
	}
	printf("\n");
#endif

#if 0
	for (int j = 0; j < ncols_; ++j)
		printf("Sol%d %e\n", j, solution_[j]);
#endif

	/** update trust region FIRST, because this does not change problem. */
	if (par_->TssDdEnableTrustRegion_)
	{
		if (newobj >= dual_bound + 1.0e-4 * (primalBound - dual_bound))
		{
			if (par_->logLevel_ > 1)
				printf("-> %s STEP: dual objective %e", isSolved_ ? "SERIOUS" : "INITIAL", newobj);

			/**
			 * 03/27/2015 - According to the paper by Linderoth and Wright, we should add cuts here.
			 */
			//if (!isSolved_)
			{
				/** add cuts and re-optimize */
				nCutsAdded = addCuts(objvals, solution);

				/** reset minor cut counter */
				ncuts_minor_ = nCutsAdded;

				/** mark cuts not to be deleted */
				for (int i = cuts_->sizeCuts() - nCutsAdded; i < cuts_->sizeCuts(); ++i)
					possiblyDelete_[i] = false;
			}

			if (isSolved_)
			{
				/** update proximal point */
#ifndef DEBUG_TRUST_REGION
#if 0
				printf("master solution:");
				for (int j = 0, jj = 0; j < ncols_; ++j)
				{
					if (fabs(solution_[j]) > 1.e-10)
					{
						if (jj > 0 && jj % 5 == 0) printf("\n");
						printf("%4d [%e] ", j, solution_[j]);
						jj++;
					}
				}
#endif
				CoinCopyN(solution_ + nCutsPerIter_, nscen_ * ncols_first_, prox_);
				if (par_->logLevel_ > 2)
					printf(", updated proximal point");
#endif

				/** possibly delete cuts */
				possiblyDeleteCuts(newobj);

#ifndef DEBUG_TRUST_REGION
				/** is solution boundary? */
				if (isSolutionBoundary() &&
					/*primalBound - dual_bound < 0 ||*/
					newobj >= dual_bound + 0.5 * (primalBound - dual_bound))
				{
					/** increase trust region */
					rho_ = CoinMin(2. * rho_, 1.0e+4);
					if (par_->logLevel_ > 2)
						printf(", increased trust region size %e", rho_);
				}
#endif

				/** set trust region */
				//printf("rho_ %e\n", rho_);
				//PRINT_ARRAY(ncols_, prox_);
				setTrustRegion(rho_, prox_);
			}

			/** update dual bound */
			dual_bound = newobj;
			trcnt_ = 0;
		}
		else
		{
#if 0
			double valid = 0.0;
			for (int s = 0; s < nscen_; ++s)
			{
				valid += objvals[s];
				for (int j = 0; j < ncols_first_; ++j)
					valid += solution[s][j] * (prox_[s * ncols_first_ + j] - solution_[nCutsPerIter_ + s * ncols_first_ + j]);
			}
			valid -= dual_bound;
			printf("%e should be positive!\n", valid);
#endif
			/** add cuts and re-optimize */
			nCutsAdded = addCuts(objvals, solution);

			/** increase minor cut counter */
			ncuts_minor_ += nCutsAdded;

			/** null step */
			if (par_->logLevel_ > 1)
				printf("-> null step: dual objective %e", newobj);

#ifndef DEBUG_TRUST_REGION
			if (primalBound < dual_bound)
			{
				/** increase trust region */
				rho_ = CoinMin(2. * rho_, 1.0e+4);
				if (par_->logLevel_ > 2)
					printf(", increased trust region size %e", rho_);
				/** set trust region */
				setTrustRegion(rho_, prox_);
			}
			else
			{
				/** The following rule is from Linderoth and Wright (2003) */
				double rho = CoinMin(1.0, rho_) * (dual_bound - newobj) / (primalBound - dual_bound);
				if (rho > 0) trcnt_++;
				if (rho >= 3 || (trcnt_ >= 3 && fabs(rho - 2.) < 1.0))
				{
					/** decrease trust region */
					rho_ *= 1.0 / CoinMin(rho, 4.);
					trcnt_ = 0;
					if (par_->logLevel_ > 2)
						printf(", decreased trust region size %e", rho_);

					/** set trust region */
					setTrustRegion(rho_, prox_);
				}
			}
#endif
		}
	}
	else
	{
		/** add cuts and re-optimize */
		nCutsAdded = addCuts(objvals, solution);

		if (par_->logLevel_ > 2)
			printf("  dual objective %e", newobj);
		if (newobj >= dual_bound + 1.0e-4 * (primalBound - dual_bound))
			/** update dual bound */
			dual_bound = newobj;
	}

	if (par_->logLevel_ > 2)
		printf(", master has %d rows and %d cols after adding %d cuts.\n",
				si_->getNumRows(), si_->getNumCols(), nCutsAdded);

#if ENABLE_UPPER_BOUNDING == 1
	/** update primal bound constraint */
	if (par_->TssDdFeasRecoveryInt_ >= 0 && si_->getColUpper()[ncols_ - 1] - primal_bound > 1.0e-10)
	{
		si_->setColUpper(ncols_ - 1, primal_bound);
		//si_->setColBounds(ncols_ - 1, primal_bound, primal_bound);
		if (par_->logLevel_ > 1)
			printf(", updated bounding cuts with rhs %e", primal_bound);
	}
#endif
	OoqpEps * ooqp = dynamic_cast<OoqpEps*>(si_);
	if (ooqp)
	{
		if (ooqp->hasOoqpStatus_ && isSolved_)
		{
			double epsilon = (si_->getPrimalBound() - newobj + ooqp->getDualityGap()) / (1. + fabs(si_->getPrimalBound()));
//				printf("epsilon %e = primalBound %e newobj %e dualityGap %e\n",
//						epsilon, si_->getPrimalBound(), newobj, ooqp->getDualityGap());
			if (epsilon > 1.) epsilon = 1.;
			ooqp->setOoqpStatus(epsilon, -primal_bound, -dual_bound);
		}
	}

	if (par_->logLevel_ > 1)
		printf("\n");

	/** set primal and dual bounds */
	primal_bounds_.push_back(primal_bound);
	dual_bounds_.push_back(dual_bound);

	return STO_RTN_OK;
}

/** solve problem */
STO_RTN_CODE TssDdPrimalMaster::solve()
{
	double stime = CoinGetTimeOfDay();

	/** solve */
	si_->solve();

	/** number of simplex iterations */
	itncnt_.push_back(si_->getIterationCount());
	itncnt_accum_ += si_->getIterationCount();

	/** mark as solved */
	isSolved_ = true;

	double etime = CoinGetTimeOfDay();

	/** solution time */
	wtime_solve_.push_back(etime - stime);
	wtime_accum_ += etime - stime;

	/** store objective value */
	if (si_->getStatus() == STO_STAT_OPTIMAL)
	{
		/** get objective */
		objval_ = si_->getPrimalBound();
		/** copy solution */
		CoinCopyN(si_->getSolution(), si_->getNumCols(), solution_);
		objvals_.push_back(objval_);
	}
	else
	{
		printf("Warning: not optimal!\n");
		DSPdebug(si_->writeMps("debug.mps"));
		objvals_.push_back(COIN_DBL_MAX);
	}

	return STO_RTN_OK;
}

/** is solution trust region boundary? */
bool TssDdPrimalMaster::isSolutionBoundary(double eps)
{
	double maxdiff = 0.0;
	double mindiff = COIN_DBL_MAX;
	const double * sol = solution_;
	const double * clbd = si_->getColLower();
	const double * cubd = si_->getColUpper();
	int ncols = si_->getNumCols();
#if ENABLE_UPPER_BOUNDING == 1
	if (par_->TssDdFeasRecoveryInt_ >= 0)
		ncols--;
#endif

	for (int j = nCutsPerIter_; j < ncols; ++j)
	{
		double diff = CoinMin(cubd[j] - sol[j], sol[j] - clbd[j]);
		if (diff > maxdiff)
			maxdiff = diff;
		if (diff < mindiff)
			mindiff = diff;
	}
	//printf("mindiff %e maxdiff %e\n", mindiff, maxdiff);
	//return fabs(maxdiff) < eps;
	return fabs(mindiff) < eps;
}

/** aging cuts */
STO_RTN_CODE TssDdPrimalMaster::agingCuts()
{
	/** exit if no cut */
	if (si_->getNumRows() == nrows_)
		return STO_RTN_OK;

	SolverInterfaceClp * clp = dynamic_cast<SolverInterfaceClp*>(si_);
	SolverInterfaceOoqp * ooqp = dynamic_cast<SolverInterfaceOoqp*>(si_);
	if (clp)
	{
		const double * pi = clp->getOSI()->getRowPrice() + nrows_;

		/** mark cuts that should not be deleted */
		for (int i = 0, i2 = 0; i < cuts_->sizeCuts(); ++i)
		{
			/** do not consider inactive cuts */
			if (cuts_age_[i] < 0) continue;

			/** aging cuts with (almost) zero Lagrangian multiplier values */
			if (fabs(pi[i2]) < 1.0e-10)
				cuts_age_[i]++;

			i2++;
		}
	}
	else if (ooqp)
	{
		const double * lambda = ooqp->lambda();
		const double * pi = ooqp->pi();

		/** mark cuts that should not be deleted */
		int numPis = ooqp->getNumPis();
		int numLambdas = ooqp->getNumLambdas();
		for (int i = cuts_->sizeCuts() - 1; i >= 0; --i)
		{
			/** do not consider inactive cuts */
			if (cuts_age_[i] < 0) continue;
			OsiRowCut * rc = cuts_->rowCutPtr(i);
			assert(rc);
			if (rc->sense() == 'G')
			{
				assert(numLambdas > 0);
				numLambdas--;
				if (fabs(lambda[numLambdas]) < 1.0e-10)
					cuts_age_[i]++;
			}
			else if (rc->sense() == 'L')
			{
				assert(numPis > 0);
				numPis--;
				if (fabs(pi[numPis]) < 1.0e-10)
					cuts_age_[i]++;
			}
		}
	}

	return STO_RTN_OK;
}

/** add cuts */
int TssDdPrimalMaster::addCuts(
		double * objvals,     /**< objective values of subproblems */
		double ** solution,   /**< subproblem solutions */
		bool      possiblyDel /**< possibly delete cuts*/)
{
#define FREE_MEMORY \
	FREE_2D_ARRAY_PTR(nCutsPerIter_, aggvec); \
	FREE_ARRAY_PTR(aggrhs);

	OsiCuts cuts;
	double ** aggvec = NULL;
	double *  aggrhs = NULL;
	CoinPackedVector cutvec; /**< cut body */
	double cutrhs;           /**< cut RHS */
#if ENABLE_UPPER_BOUNDING == 2
	CoinPackedVector cutvec2; /**< cut body */
	double cutrhs2 = 0.0;     /**< cut RHS */
#endif

#if 1
	/** allocate memory for dense cut */
	aggvec = new double * [nCutsPerIter_];
	aggrhs = new double [nCutsPerIter_];
	for (int i = 0; i < nCutsPerIter_; ++i)
	{
		aggvec[i] = new double [ncols_];
		CoinZeroN(aggvec[i], ncols_);
		aggrhs[i] = 0.;
	}

	/** add row cuts by looping over each scenario */
	for (int s = 0; s < nscen_; ++s)
	{
		assert(fabs(objvals[s]) < 1.0e+20);

		/** cut index */
		int cutidx = s % nCutsPerIter_;

		/** construct cut */
		aggrhs[cutidx] += objvals[s];
		for (int j = 0; j < ncols_first_; ++j)
		{
			if (fabs(solution[s][j]) > 1E-10)
			{
				aggvec[cutidx][nCutsPerIter_ + s * ncols_first_ + j] -= solution[s][j];
				if (isSolved_)
					aggrhs[cutidx] -= solution_[nCutsPerIter_ + s * ncols_first_ + j] * solution[s][j];
			}
		}
	}

#if 0
	printf("\n");
	for (int i = 0; i < cuts_->sizeCuts(); ++i)
	{
		OsiRowCut *    oldcut = cuts_->rowCutPtr(i);
		double         tmpsum = oldcut->rhs();
		int            oldnum = oldcut->row().getNumElements();
		const int *    oldind = oldcut->row().getIndices();
		const double * oldval = oldcut->row().getElements();
		for (int j = 0; j < oldnum; ++j)
		{
			if (oldind[j] < nCutsPerIter_) continue;
			tmpsum -= oldval[j] * prox_[oldind[j]-nCutsPerIter_];
		}
		printf("cut(%d): theta <= %e at proximal point\n", i, tmpsum);
	}
#endif

	for (int s = 0; s < nCutsPerIter_; ++s)
	{
		aggvec[s][s] = 1.0;

		/** construct cut */
		cutvec.clear();

#if 0
		double coefsum = 0.0;
		double coefmax = -COIN_DBL_MAX;
		double coefmin = COIN_DBL_MAX;
#endif

		/** set it as sparse */
		for (int j = 0; j < ncols_; ++j)
		{
#if 0
			coefsum += aggvec[s][j];
			if (coefmax < aggvec[s][j]) coefmax = aggvec[s][j];
			if (coefmin > aggvec[s][j]) coefmin = aggvec[s][j];
#endif
			if (fabs(aggvec[s][j]) > 1E-10)
				cutvec.insert(j, aggvec[s][j]);
		}

		/** cut rhs */
		cutrhs = aggrhs[s];
		if (fabs(cutrhs) < 1E-10)
			cutrhs = 0.0;

		OsiRowCut rc;
		rc.setRow(cutvec);
		rc.setLb(-COIN_DBL_MAX);
		rc.setUb(cutrhs);
		if (isSolved_)
			rc.setEffectiveness(rc.violated(solution_) / cutvec.twoNorm());
		else
			rc.setEffectiveness(1.0);
		//rc.print();
#if 0
		double tmpsum = cutrhs;
		for (int j = nCutsPerIter_; j < ncols_; ++j)
			tmpsum -= aggvec[s][j] * prox_[j-nCutsPerIter_];
		printf("new cut: theta(%d) <= %e at proximal point\n", s, tmpsum);
		printf("Cut: rhs %e nzcnt %d sum %e avg %e min %e max %e\n",
				rc.ub(), cutvec.getNumElements(), coefsum, coefsum / ncols_, coefmin, coefmax);
		printf("violation: %e\n", rc.violated(solution_));
#endif

		if (rc.effectiveness() > 1.0E-6)
		{
			/** number of cuts before adding cut */
			int nCutsBefore = cuts_->sizeCuts();

			/** add cut if not duplicate */
			cuts_->insertIfNotDuplicate(rc);

			if (nCutsBefore < cuts_->sizeCuts())
			{
				/** insertIfNotDuplicate does not set effectiveness */
				cuts_->rowCutPtr(nCutsBefore)->setEffectiveness(rc.effectiveness());
				cuts_age_.push_back(0);
				possiblyDelete_.push_back(possiblyDel);
				master_objvals_.push_back(si_->getPrimalBound());
				cuts.insert(rc);
			}
		}
	}
#else
	/** add row cuts by looping over each scenario */
	for (int s = 0; s < nscen_; ++s)
	{
		assert(fabs(objvals[s]) < 1.0e+20);
		/** construct cut */
		cutvec.clear();
		cutvec.insert(s, 1.0);
		cutrhs = objvals[s];
#if ENABLE_UPPER_BOUNDING == 2
		if (par_->TssDdFeasRecoveryInt_ >= 0)
			cutrhs2 += objvals[s];
#endif
		for (int j = 0; j < ncols_first_; ++j)
		{
			if (fabs(solution[s][j]) > 1E-10)
			{
				cutvec.insert(nscen_ + s * ncols_first_ + j, -solution[s][j]);
#if ENABLE_UPPER_BOUNDING == 2
				if (par_->TssDdFeasRecoveryInt_ >= 0)
					cutvec2.insert(nscen_ + s * ncols_first_ + j, -solution[s][j]);
#endif
				if (isSolved_)
				{
					cutrhs -= solution_[nscen_ + s * ncols_first_ + j] * solution[s][j];
#if ENABLE_UPPER_BOUNDING == 2
					if (par_->TssDdFeasRecoveryInt_ >= 0)
						cutrhs2 -= solution_[nscen_ + s * ncols_first_ + j] * solution[s][j];
#endif
				}
			}
		}

		OsiRowCut rc;
		rc.setRow(cutvec);
		rc.setLb(-COIN_DBL_MAX);
		rc.setUb(cutrhs);
		if (isSolved_)
			rc.setEffectiveness(rc.violated(solution_));
		else
			rc.setEffectiveness(1.0);
		//rc.print();

		if (rc.effectiveness() > 1.0E-10)
		{
			/** number of cuts before adding cut */
			int nCutsBefore = cuts_->sizeCuts();

			/** add cut if not duplicate */
			cuts_->insertIfNotDuplicate(rc);

			if (nCutsBefore < cuts_->sizeCuts())
			{
				/** insertIfNotDuplicate does not set effectiveness */
				cuts_->rowCutPtr(nCutsBefore)->setEffectiveness(rc.effectiveness());
				cuts_age_.push_back(0);
				possiblyDelete_.push_back(possiblyDel);
				master_objvals_.push_back(si_->getPrimalBound());
				cuts.insert(rc);
			}
		}
	}
#endif

#if ENABLE_UPPER_BOUNDING == 2
	/** add bounding cuts */
	if (par_->TssDdFeasRecoveryInt_ >= 0)
	{
		cutvec2.insert(ncols_ - 1, 1.0);
		OsiRowCut rc;
		rc.setRow(cutvec2);
		rc.setLb(cutrhs2);
		rc.setUb(COIN_DBL_MAX);
//		if (isSolved_)
//			rc.setEffectiveness(rc.violated(solution_));
//		else
			rc.setEffectiveness(1.0);
		//rc.print();
		//printf("lhs %e >= rhs %e\n", rc.row().dotProduct(si_->getSolution()), cutrhs2);
		if (rc.effectiveness() > 1.0E-10)
			ubCuts_->insert(rc);
	}
#endif

	FREE_MEMORY

	/** exit if no cut is found */
	if (cuts.sizeCuts() == 0)
		return 0;

	/** apply cuts */
	int nCutsAdded = si_->addCuts(cuts);
	if (nCutsAdded == 0)
	{
		/** recruit back some cuts if no cut is generated */
		recruiteCuts();
	}

	return nCutsAdded;
#undef FREE_MEMORY
}

/** possibly delete cuts */
STO_RTN_CODE TssDdPrimalMaster::possiblyDeleteCuts(
		double sub_objval /**< subproblem objective value */)
{
	OsiCuts cuts;
	int ncuts = si_->getNumRows() - nrows_;

	if (ncuts == 0)
		return STO_RTN_OK;

	SolverInterfaceOsi * osi = dynamic_cast<SolverInterfaceOsi*>(si_);
	SolverInterfaceOoqp * ooqp = dynamic_cast<SolverInterfaceOoqp*>(si_);
	if (osi)
	{
		const double * pi = osi->getOSI()->getRowPrice() + nrows_;

		/** mark cuts that should not be deleted */
		for (int i = 0, i2 = 0; i < cuts_->sizeCuts(); ++i)
		{
			/** do not consider inactive cuts */
			if (cuts_age_[i] < 0) continue;

			/** consider only old cuts */
			if (cuts_age_[i] < 100)
			{
				possiblyDelete_[i] = false;
				continue;
			}

			/** aging cuts with (almost) zero Lagrangian multiplier values */
			if (fabs(pi[i2]) < 1.0e-10)
				possiblyDelete_[i] = false;
			/** do not delete cuts generated at minor iterations such that the following condition holds. */
			else if (i >= cuts_->sizeCuts() - ncuts_minor_ &&
					(si_->getPrimalBound() - sub_objval) > cutdel_param_ * (master_objvals_[i] - sub_objval))
				possiblyDelete_[i] = false;

			i2++;
		}
	}
	else if (ooqp)
	{
		const double * lambda = ooqp->lambda();
		const double * pi = ooqp->pi();

		/** mark cuts that should not be deleted */
		int numPis = ooqp->getNumPis();
		int numLambdas = ooqp->getNumLambdas();
		for (int i = cuts_->sizeCuts() - 1; i >= 0; --i)
		{
			/** do not consider inactive cuts */
			if (cuts_age_[i] < 0) continue;

			/** consider only old cuts */
			if (cuts_age_[i] < 100)
			{
				possiblyDelete_[i] = false;
				continue;
			}

			OsiRowCut * rc = cuts_->rowCutPtr(i);
			assert(rc);
			if (rc->sense() == 'G')
			{
				numLambdas--;
				if (fabs(lambda[numLambdas]) < 1.0e-10)
					possiblyDelete_[i] = false;
				/** do not delete cuts generated at minor iterations such that the following condition holds. */
				else if (i >= cuts_->sizeCuts() - ncuts_minor_ &&
						(si_->getPrimalBound() - sub_objval) > cutdel_param_ * (master_objvals_[i] - sub_objval))
					possiblyDelete_[i] = false;
			}
			else if (rc->sense() == 'L')
			{
				numPis--;
				if (fabs(pi[numPis]) < 1.0e-10)
					possiblyDelete_[i] = false;
				/** do not delete cuts generated at minor iterations such that the following condition holds. */
				else if (i >= cuts_->sizeCuts() - ncuts_minor_ &&
						(si_->getPrimalBound() - sub_objval) > cutdel_param_ * (master_objvals_[i] - sub_objval))
					possiblyDelete_[i] = false;
			}
		}
	}

	/** get basis information */
	CoinWarmStartBasis * ws = NULL;
	if (osi) ws = dynamic_cast<CoinWarmStartBasis*>(osi->getWarmStart()->clone());

	vector<char> aStat; /**< status of artificial variables */

	/** mark as deleted; and construct temporary cut pool to be added */
	for (int i = 0, i2 = nrows_; i < cuts_->sizeCuts(); ++i)
	{
		/** do not consider inactive cuts */
		if (cuts_age_[i] < 0) continue;

		if (possiblyDelete_[i])
			cuts_age_[i] = -1;
		else
		{
			OsiRowCut * rc = cuts_->rowCutPtr(i);
			if (rc)
			{
				rc->setEffectiveness(1.0);
				cuts.insert(*rc);
				if (osi)
					aStat.push_back(ws->getArtifStatus(i2));
			}
		}

		i2++;
	}

	/** number of cuts to delete */
	int nCutsToDelete = ncuts - cuts.sizeCuts();
	//cout << ", found " << nCutsToDelete << " cuts to delete";

	/** exit if no cut to delete */
	if (nCutsToDelete == 0)
		return STO_RTN_OK;

	/** remove all cuts from solver interface */
	removeAllCuts();

	/** apply cuts */
	si_->addCuts(cuts);

	if (osi)
	{
		/** create new basis */
		CoinWarmStartBasis * basis = new CoinWarmStartBasis(
				ws->getNumStructural(), ws->getNumArtificial(),
				ws->getStructuralStatus(), &aStat[0]);

		osi->setWarmStart(basis);
	}

	return STO_RTN_OK;
}

/** recruit cuts: investigating and adding effective cuts (that were deleted) in the cut pool */
int TssDdPrimalMaster::recruiteCuts()
{
	int nRecruited = 0;
	OsiCuts cuts;

	CoinWarmStartBasis * ws = NULL;
	SolverInterfaceClp * clp = dynamic_cast<SolverInterfaceClp*>(si_);

	vector<char> aStat; /**< status of artificial variables */

	if (clp)
	{
		/** get basis information */
		clp->setWarmStart(clp->getOSI()->getWarmStart());
		ws = dynamic_cast<CoinWarmStartBasis*>(clp->getWarmStart());
	}

	//cout << ", investigating removed cuts ...";
	for (int i = 0, irow = nrows_; i < cuts_->sizeCuts(); ++i)
	{
		/** retrieve row cut */
		OsiRowCut * rc = cuts_->rowCutPtr(i);
		assert(rc);

		if (cuts_age_[i] >= 0)
		{
			/** add cut */
			cuts.insert(*rc);
			if (clp)
			{
				/** set status of artificial variable */
				aStat.push_back(ws->getArtifStatus(irow++));
			}
		}
		else
		{
			/** set effectiveness */
			rc->setEffectiveness(rc->violated(solution_));
			if (rc->effectiveness() > 1.0e-10)
			{
				nRecruited++;
				/** add cut */
				cuts.insert(*rc);
				if (clp)
				{
					/** set status of artificial variable */
					aStat.push_back(CoinWarmStartBasis::basic);
				}

				/** other cut info */
				cuts_age_[i] = 0;
				possiblyDelete_[i] = true;
			}
		}
	}
	//cout << " added back " << nRecruited << " cuts";

	if (cuts.sizeCuts() == 0)
		return 0;

	/** remove all cuts from solver interface */
	removeAllCuts();

	/** apply cuts */
	si_->addCuts(cuts);

	if (clp)
	{
		/** create new basis */
		CoinWarmStartBasis * basis = new CoinWarmStartBasis(
				ws->getNumStructural(), ws->getNumArtificial(),
				ws->getStructuralStatus(), &aStat[0]);
		clp->setWarmStart(basis);
	}

	return nRecruited;
}

/** remove all cuts from solver interface */
STO_RTN_CODE TssDdPrimalMaster::removeAllCuts()
{
	SolverInterfaceClp * clp = dynamic_cast<SolverInterfaceClp*>(si_);
	SolverInterfaceOoqp * ooqp = dynamic_cast<SolverInterfaceOoqp*>(si_);
	if (clp)
	{
		int ncuts = si_->getNumRows() - nrows_;

		/** row indices to delete */
		int * rowIndices = new int [ncuts];
		CoinIotaN(rowIndices, ncuts, nrows_);

		/** delete */
		clp->getOSI()->deleteRows(ncuts, rowIndices);

		/** free memory */
		FREE_ARRAY_PTR(rowIndices);
	}
	else if (ooqp)
	{
		for (int i = 0; i < ooqp->cuts_.sizeCuts(); ++i)
			delete ooqp->cuts_.rowCutPtr(i);
		ooqp->cuts_.dumpCuts();
	}

	return STO_RTN_OK;
}

