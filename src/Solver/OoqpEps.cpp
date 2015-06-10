/*
 * OoqpEps.cpp
 *
 *  Created on: Jan 28, 2015
 *      Author: kibaekkim
 */

/** DSP */
#include "Solver/OoqpEps.h"
#include "Solver/OoqpStatus.h"
#include "Utility/StoMacros.h"

/** set my OOQP status */
void OoqpEps::setOoqpStatus(double epsilon, double lowerBound, double upperBound)
{
	releaseOOQP();
	epsilon_       = epsilon;
	lowerBound_    = lowerBound;
	upperBound_    = upperBound;
}

/** solve */
void OoqpEps::solve()
{
	if (released_)
		gutsOfLoadProblem();

	//prob_->print();
	if (print_level_ > 0) solver_->monitorSelf();
	int status = solver_->solve(prob_, vars_, resid_);
	switch(status)
	{
	case 0:
	{
		status_ = STO_STAT_OPTIMAL;
		objval_ = prob_->objectiveValue(vars_);

		/** get variable values */
		vars_->x->copyIntoArray(x_);
		vars_->y->copyIntoArray(y_);
		vars_->lambda->copyIntoArray(lambda_);
		vars_->pi->copyIntoArray(pi_);
		vars_->gamma->copyIntoArray(gamma_);
		vars_->phi->copyIntoArray(phi_);

		nLambdas_ =  vars_->lambda->length();
		nPis_ = vars_->pi->length();
		dualityGap_ = resid_->dualityGap();

		/** number of iterations */
		nIters_ = solver_->iter;

		/** suboptimal? */
		if (vars_->mu() <= solver_->getMuTol() &&
			resid_->residualNorm() <= solver_->getArTol() * prob_->datanorm())
			suboptimal_ = false;
		else
			suboptimal_ = true;
		//printf("optimal? %s\n", suboptimal_ ? "no" : "yes");

		break;
	}
	case 1:
		status_ = STO_STAT_STOPPED_UNKNOWN;
		break;
	case 2:
		status_ = STO_STAT_STOPPED_ITER;
		break;
	case 3:
		status_ = STO_STAT_PRIM_INFEASIBLE;
		break;
	case 4:
	default:
		status_ = STO_STAT_UNKNOWN;
		break;
	}
}

/** core part for load problem */
void OoqpEps::gutsOfLoadProblem()
{
#define FREE_MEMORY \
		FREE_ARRAY_PTR(c); \
		FREE_ARRAY_PTR(xlow); \
		FREE_ARRAY_PTR(ixlow); \
		FREE_ARRAY_PTR(xupp); \
		FREE_ARRAY_PTR(ixupp); \
		FREE_ARRAY_PTR(irowA); \
		FREE_ARRAY_PTR(jcolA); \
		FREE_ARRAY_PTR(dA); \
		FREE_ARRAY_PTR(b); \
		FREE_ARRAY_PTR(irowC); \
		FREE_ARRAY_PTR(jcolC); \
		FREE_ARRAY_PTR(dC); \
		FREE_ARRAY_PTR(clow); \
		FREE_ARRAY_PTR(iclow); \
		FREE_ARRAY_PTR(cupp); \
		FREE_ARRAY_PTR(icupp);

	int nx = 0;
	int nnzA = 0;
	int nnzC = 0;
	double * c     = NULL;
	double * xlow  = NULL;
	char *   ixlow = NULL;
	double * xupp  = NULL;
	char *   ixupp = NULL;
	int *    irowA = NULL;
	int *    jcolA = NULL;
	double * dA    = NULL;
	double * b     = NULL;
	int *    irowC = NULL;
	int *    jcolC = NULL;
	double * dC    = NULL;
	double * clow  = NULL;
	char *   iclow = NULL;
	double * cupp  = NULL;
	char *   icupp = NULL;

	BGN_TRY_CATCH

	/** number of variables */
	nx = ncols_;

	/** count number of equality/inequality constraints */
	/** count number of nonzero elements in A and C */
	my_ = 0;
	mz_ = 0;
	for (int i = 0; i < mat_->getNumRows(); ++i)
	{
		if (fabs(rlbd_[i] - rubd_[i]) < 1.0e-10)
		{
			my_++;
			nnzA += mat_->getVectorSize(i);
		}
		else
		{
			mz_++;
			nnzC += mat_->getVectorSize(i);
		}
	}
	for (int i = 0; i < cuts_.sizeCuts(); ++i)
	{
		OsiRowCut * rc = cuts_.rowCutPtr(i);
		if (!rc) continue;
		if (rc->lb() == rc->ub())
		{
			my_++;
			nnzA += rc->row().getNumElements();
		}
		else
		{
			mz_++;
			nnzC += rc->row().getNumElements();
		}
	}

	/** linear objective */
	c = new double [nx];
	for (int j = 0; j < nx; ++j)
		c[j] = obj_[j] * sense_;

	/** column bounds */
	xlow  = new double [nx];
	ixlow = new char [nx];
	xupp  = new double [nx];
	ixupp = new char [nx];
	for (int j = 0; j < nx; ++j)
	{
		if (clbd_[j] > -COIN_DBL_MAX)
		{
			xlow[j] = clbd_[j];
			ixlow[j] = 1;
		}
		else
		{
			xlow[j] = 0;
			ixlow[j] = 0;
		}
		if (cubd_[j] < COIN_DBL_MAX)
		{
			xupp[j] = cubd_[j];
			ixupp[j] = 1;
		}
		else
		{
			xupp[j] = 0;
			ixupp[j] = 0;
		}
	}

	/** for constraints */
	irowA = new int [nnzA];
	jcolA = new int [nnzA];
	dA    = new double [nnzA];
	b     = new double [my_];
	irowC = new int [nnzC];
	jcolC = new int [nnzC];
	dC    = new double [nnzC];
	clow  = new double [mz_];
	iclow = new char [mz_];
	cupp  = new double [mz_];
	icupp = new char [mz_];
	int posA = 0;
	int posC = 0;
	int posy = 0;
	int posz = 0;
	for (int i = 0; i < mat_->getNumRows(); ++i)
	{
		if (fabs(rlbd_[i] - rubd_[i]) < 1.0e-10)
		{
			for (int j = 0; j < mat_->getVectorSize(i); ++j)
			{
				irowA[posA] = posy;
				jcolA[posA] = mat_->getIndices()[mat_->getVectorFirst(i) + j];
				dA[posA] = mat_->getElements()[mat_->getVectorFirst(i) + j];
				posA++;
			}
			b[posy] = rubd_[i];
			posy++;
		}
		else
		{
			for (int j = 0; j < mat_->getVectorSize(i); ++j)
			{
				irowC[posC] = posz;
				jcolC[posC] = mat_->getIndices()[mat_->getVectorFirst(i) + j];
				dC[posC] = mat_->getElements()[mat_->getVectorFirst(i) + j];
				posC++;
			}
			if (rlbd_[i] > -COIN_DBL_MAX)
			{
				clow[posz]  = rlbd_[i];
				iclow[posz] = 1;
			}
			else
			{
				clow[posz]  = 0;
				iclow[posz] = 0;
			}
			if (rubd_[i] < COIN_DBL_MAX)
			{
				cupp[posz]  = rubd_[i];
				icupp[posz] = 1;
			}
			else
			{
				cupp[posz]  = 0;
				icupp[posz] = 0;
			}
			posz++;
		}
	}
	for (int i = 0; i < cuts_.sizeCuts(); ++i)
	{
		OsiRowCut * rc = cuts_.rowCutPtr(i);
		if (!rc) continue;
		const CoinPackedVector row = rc->row();
		if (rc->lb() == rc->ub())
		{
			for (int j = 0; j < row.getNumElements(); ++j)
			{
				irowA[posA] = posy;
				jcolA[posA] = row.getIndices()[j];
				dA[posA] = row.getElements()[j];
				posA++;
			}
			b[posy] = rc->ub();
			posy++;
		}
		else
		{
			for (int j = 0; j < row.getNumElements(); ++j)
			{
				irowC[posC] = posz;
				jcolC[posC] = row.getIndices()[j];
				dC[posC] = row.getElements()[j];
				posC++;
			}
			if (rc->lb() > -COIN_DBL_MAX)
			{
				clow[posz]  = rc->lb();
				iclow[posz] = 1;
			}
			else
			{
				clow[posz]  = 0;
				iclow[posz] = 0;
			}
			if (rc->ub() < COIN_DBL_MAX)
			{
				cupp[posz]  = rc->ub();
				icupp[posz] = 1;
			}
			else
			{
				cupp[posz]  = 0;
				icupp[posz] = 0;
			}
			posz++;
		}
	}

	FREE_ARRAY_PTR(x_);
	FREE_ARRAY_PTR(y_);
	FREE_ARRAY_PTR(lambda_);
	FREE_ARRAY_PTR(pi_);
	FREE_ARRAY_PTR(gamma_);
	FREE_ARRAY_PTR(phi_);
	x_      = new double [nx];
	y_      = new double [my_];
	lambda_ = new double [mz_];
	pi_     = new double [mz_];
	gamma_  = new double [nx];
	phi_    = new double [nx];

#ifdef DSP_HAS_MA57
	qp_ = new QpGenSparseMa57(nx, my_, mz_, 0, nnzA, nnzC);
#else
	qp_ = new QpGenSparseMa27(nx, my_, mz_, 0, nnzA, nnzC);
#endif

	/**
	 * Data storage for the regularized master problem:
	 *   Q is the lower triangular matrix of the quadratic term of the objective function.
	 *   A is the equality constraint matrix.
	 *   C is the inequality constraint matrix.
	 */
	prob_ = (QpGenData*)qp_->copyDataFromSparseTriple(
			c,     NULL,  0,     NULL,  NULL,
			xlow,  ixlow, xupp,  ixupp,
			irowA, nnzA,  jcolA, dA,    b,
			irowC, nnzC,  jcolC, dC,
			clow,  iclow, cupp,  icupp);
	//prob_->print();

	/** declare variables */
	vars_ = (QpGenVars*)qp_->makeVariables(prob_);

	/** declare residuals */
	resid_ = (QpGenResiduals*)qp_->makeResiduals(prob_);

	/** create solver */
	//solver_ = new GondzioSolver(qp_, prob_);
	solver_ = new MehrotraSolver(qp_, prob_);

	if (hasOoqpStatus_)
	{
		OoqpStatus * mystat = new OoqpStatus(epsilon_, lowerBound_, upperBound_);
		solver_->useStatus(mystat);
	}

	released_ = false;

	END_TRY_CATCH(FREE_MEMORY)

	FREE_MEMORY

#undef FREE_MEMORY
}

