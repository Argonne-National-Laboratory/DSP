/*
 * SolverInterfaceOoqp.cpp
 *
 *  Created on: Jan 7, 2015
 *      Author: kibaekkim
 */

//#define DSP_DEBUG

/** DSP */
#include <Utility/DspMacros.h>
#include <Utility/DspMpi.h>
#include "SolverInterface/SolverInterfaceOoqp.h"
#include "GondzioSolver.h"
#include "MehrotraSolver.h"
//#ifdef DSP_HAS_MA57
//#include "QpGenSparseMa57.h"
//#else
//#endif
#include "QpGenSparseMa27.h"
#include "QpGenSparseMa57.h"
#include "QpGenDense.h"

/** copy constructor */
SolverInterfaceOoqp::SolverInterfaceOoqp(SolverInterfaceOoqp * si) :
	SolverInterface(si->par_)
{
	if (si->dynCols_)
		dynCols_ = new dynColumns(si->dynCols_);

	sense_ = si->getObjSense();
	status_ = si->getStatus();
	print_level_ = si->getPrintLevel();
	released_ = si->released();

	my_ = si->getNumEqConss();
	mz_ = si->getNumIneqConss();
	nrows_    = my_ + mz_;
	ncols_    = si->getNumCols();
	nLambdas_ = si->getNumLambdas();
	nPis_     = si->getNumPis();
	objval_   = si->getPrimalBound() * sense_;
	nIters_   = si->getIterationCount();

	/** copy matrix */
	mat_ = new CoinPackedMatrix(*(si->getMatrix()));

	/** allocate memory */
	clbd_   = new double [ncols_];
	cubd_   = new double [ncols_];
	obj_    = new double [ncols_];
	rlbd_   = new double [nrows_];
	rubd_   = new double [nrows_];
	x_      = new double [ncols_];
	y_      = new double [my_];
	lambda_ = new double [nLambdas_];
	pi_     = new double [nPis_];
	gamma_  = new double [ncols_];
	phi_    = new double [ncols_];

	/** copy data */
	CoinCopyN(si->getColLower(), ncols_, clbd_);
	CoinCopyN(si->getColUpper(), ncols_, cubd_);
	CoinCopyN(si->getObjCoef(), ncols_, obj_);
	CoinCopyN(si->getRowLower(), nrows_, rlbd_);
	CoinCopyN(si->getRowUpper(), nrows_, rubd_);
	CoinCopyN(si->getSolution(), ncols_, x_);
	CoinCopyN(si->y(), my_, y_);
	CoinCopyN(si->lambda(), nLambdas_, lambda_);
	CoinCopyN(si->pi(), nPis_, pi_);
	CoinCopyN(si->gamma(), ncols_, gamma_);
	CoinCopyN(si->phi(), ncols_, phi_);

	releaseOOQP();
	GUTS_OF_LOAD_PROBLEM();
}

/** clone */
SolverInterface * SolverInterfaceOoqp::clone()
{
	return new SolverInterfaceOoqp(this);
}

SolverInterfaceOoqp::~SolverInterfaceOoqp()
{
	BGN_TRY_CATCH
	DSP_RTN_CHECK_THROW(finalize());
	END_TRY_CATCH(;)
}

/** finalize solver interface */
DSP_RTN_CODE SolverInterfaceOoqp::finalize()
{
	releaseOOQP();
	release();
	return DSP_RTN_OK;
}

/** release OOQP objects */
void SolverInterfaceOoqp::releaseOOQP()
{
	FREE_PTR(solver_);
	FREE_PTR(resid_);
	FREE_PTR(vars_);
	FREE_PTR(prob_);
	FREE_PTR(qp_);
	released_ = true;
}

/** release solver */
void SolverInterfaceOoqp::release()
{
	BGN_TRY_CATCH

	FREE_PTR(dynCols_);
	FREE_PTR(mat_);
	FREE_ARRAY_PTR(xlow_);
	FREE_ARRAY_PTR(ixlow_);
	FREE_ARRAY_PTR(xupp_);
	FREE_ARRAY_PTR(ixupp_);
	FREE_ARRAY_PTR(c_);
	FREE_ARRAY_PTR(clbd_);
	FREE_ARRAY_PTR(cubd_);
	FREE_ARRAY_PTR(obj_);
	FREE_ARRAY_PTR(rlbd_);
	FREE_ARRAY_PTR(rubd_);
	FREE_ARRAY_PTR(irowA_);
	FREE_ARRAY_PTR(jcolA_);
	FREE_ARRAY_PTR(dA_);
	FREE_ARRAY_PTR(b_);
	FREE_ARRAY_PTR(irowC_);
	FREE_ARRAY_PTR(jcolC_);
	FREE_ARRAY_PTR(dC_);
	FREE_ARRAY_PTR(clow_);
	FREE_ARRAY_PTR(iclow_);
	FREE_ARRAY_PTR(cupp_);
	FREE_ARRAY_PTR(icupp_);
	FREE_ARRAY_PTR(irowQ_);
	FREE_ARRAY_PTR(jcolQ_);
	FREE_ARRAY_PTR(dQ_);
	FREE_ARRAY_PTR(x_);
	FREE_ARRAY_PTR(y_);
	FREE_ARRAY_PTR(lambda_);
	FREE_ARRAY_PTR(pi_);
	FREE_ARRAY_PTR(gamma_);
	FREE_ARRAY_PTR(phi_);
	nrows_ = 0;
	ncols_ = 0;

	END_TRY_CATCH(;)
}

/** load problem */
void SolverInterfaceOoqp::loadProblem(
		OsiSolverInterface * si,
		const char * probname)
{
	CoinPackedMatrix * mat = new CoinPackedMatrix(*si->getMatrixByRow());
	loadProblem(mat,
			si->getColLower(), si->getColUpper(),
			si->getObjCoefficients(),
			NULL, si->getRowLower(), si->getRowUpper(), probname);
	FREE_PTR(mat);
}

/** load problem */
void SolverInterfaceOoqp::loadProblem(
		CoinPackedMatrix * mat,
		const double * collb,
		const double * colub,
		const double * obj,
		const char * ctype,
		const double * rowlb,
		const double * rowub,
		const char * probname)
{
	/** release everything */
	releaseOOQP();
	release();

	/** want matrix to be row-wise */
	if (mat->isColOrdered()) mat->reverseOrdering();

	/** save original data */
	mat_ = new CoinPackedMatrix(*mat);
	ncols_ = mat_->getNumCols();
	nrows_ = mat_->getNumRows();
	clbd_ = new double [ncols_];
	cubd_ = new double [ncols_];
	obj_  = new double [ncols_];
	rlbd_ = new double [nrows_];
	rubd_ = new double [nrows_];
	CoinCopyN(collb, ncols_, clbd_);
	CoinCopyN(colub, ncols_, cubd_);
	CoinCopyN(obj, ncols_, obj_);
	CoinCopyN(rowlb, nrows_, rlbd_);
	CoinCopyN(rowub, nrows_, rubd_);

	GUTS_OF_LOAD_PROBLEM();
}

void SolverInterfaceOoqp::gutsOfLoadDenseProblem()
{
	int ndynCols = 0;
	CoinPackedMatrix * cols = NULL; /**< do not free */

	BGN_TRY_CATCH

	/** retrieve dynamic column information */
	if (dynCols_)
	{
		ndynCols = dynCols_->size();
		cols = dynCols_->cols_;
	}

	/** number of variables */
	nx_ = ncols_ + ndynCols;

	/** count number of equality/inequality constraints */
	/** count number of nonzero elements in A and C */
	my_ = 0;
	mz_ = 0;
	for (int i = 0; i < mat_->getNumRows(); ++i)
	{
		if (fabs(rlbd_[i] - rubd_[i]) < 1.0e-10)
			my_++;
		else
			mz_++;
	}
	for (int i = 0; i < cuts_.sizeCuts(); ++i)
	{
		OsiRowCut * rc = cuts_.rowCutPtr(i);
		if (!rc) continue;
		if (rc->lb() == rc->ub())
			my_++;
		else
			mz_++;
	}

	/** linear objective */
	FREE_ARRAY_PTR(c_);
	c_ = new double [nx_];
	for (int j = 0; j < ncols_; ++j)
		c_[j] = obj_[j] * sense_;
	/** for dynamic columns */
	for (int j = ncols_; j < nx_; ++j)
		c_[j] = dynCols_->objs_[j-ncols_] * sense_;
	//DspMessage::printArray(nx_, c_);

	/** column bounds */
	FREE_ARRAY_PTR(xlow_);
	FREE_ARRAY_PTR(ixlow_);
	FREE_ARRAY_PTR(xupp_);
	FREE_ARRAY_PTR(ixupp_);
	xlow_  = new double [nx_];
	ixlow_ = new char [nx_];
	xupp_  = new double [nx_];
	ixupp_ = new char [nx_];
	for (int j = 0; j < ncols_; ++j)
	{
		if (clbd_[j] > -COIN_DBL_MAX)
		{
			xlow_[j] = clbd_[j];
			ixlow_[j] = 1;
		}
		else
		{
			xlow_[j] = 0;
			ixlow_[j] = 0;
		}
		if (cubd_[j] < COIN_DBL_MAX)
		{
			xupp_[j] = cubd_[j];
			ixupp_[j] = 1;
		}
		else
		{
			xupp_[j] = 0;
			ixupp_[j] = 0;
		}
	}

	/** for dynamic columns */
	for (int j = ncols_; j < nx_; ++j)
	{
		if (dynCols_->lbs_[j-ncols_] > -COIN_DBL_MAX)
		{
			xlow_[j] = dynCols_->lbs_[j-ncols_];
			ixlow_[j] = 1;
		}
		else
		{
			xlow_[j] = 0;
			ixlow_[j] = 0;
		}
		if (dynCols_->ubs_[j-ncols_] < COIN_DBL_MAX)
		{
			xupp_[j] = dynCols_->ubs_[j-ncols_];
			ixupp_[j] = 1;
		}
		else
		{
			xupp_[j] = 0;
			ixupp_[j] = 0;
		}
	}

	/** for constraints */
	FREE_ARRAY_PTR(dA_);
	FREE_ARRAY_PTR(b_);
	FREE_ARRAY_PTR(dC_);
	FREE_ARRAY_PTR(clow_);
	FREE_ARRAY_PTR(iclow_);
	FREE_ARRAY_PTR(cupp_);
	FREE_ARRAY_PTR(icupp_);
	dA_    = new double [nx_*my_];
	b_     = new double [my_];
	dC_    = new double [nx_*mz_];
	clow_  = new double [mz_];
	iclow_ = new char [mz_];
	cupp_  = new double [mz_];
	icupp_ = new char [mz_];
	CoinZeroN(dA_, nx_*my_);
	CoinZeroN(dC_, nx_*mz_);
	int posA = 0, posC = 0, posy = 0, posz = 0;
	for (int i = 0; i < mat_->getNumRows(); ++i)
	{
		if (fabs(rlbd_[i] - rubd_[i]) < 1.0e-10)
		{
			for (int j = 0; j < mat_->getVectorSize(i); ++j)
			{
				int ind = posy * nx_ + mat_->getIndices()[mat_->getVectorFirst(i) + j];
				dA_[ind] = mat_->getElements()[mat_->getVectorFirst(i) + j];
				posA++;
			}
			b_[posy] = rubd_[i];
			posy++;
		}
		else
		{
			for (int j = 0; j < mat_->getVectorSize(i); ++j)
			{
				int ind = posz * nx_ + mat_->getIndices()[mat_->getVectorFirst(i) + j];
				dC_[ind] = mat_->getElements()[mat_->getVectorFirst(i) + j];
				posC++;
			}
			if (rlbd_[i] > -COIN_DBL_MAX)
			{
				clow_[posz]  = rlbd_[i];
				iclow_[posz] = 1;
			}
			else
			{
				clow_[posz]  = 0;
				iclow_[posz] = 0;
			}
			if (rubd_[i] < COIN_DBL_MAX)
			{
				cupp_[posz]  = rubd_[i];
				icupp_[posz] = 1;
			}
			else
			{
				cupp_[posz]  = 0;
				icupp_[posz] = 0;
			}
			posz++;
		}
		/** for dynamic columns */
		for (int j = 0; j < ndynCols; ++j)
		{
			for (int rowi = cols->getVectorFirst(j); rowi < cols->getVectorLast(j); ++rowi)
			{
				if (i != cols->getIndices()[rowi])
					continue;
				//DSPdebugMessage("irowA %d jcolA %d dA %e\n", i, ncols_ + j, cols->getElements()[rowi]);
				if (fabs(rlbd_[i] - rubd_[i]) < 1.0e-10)
				{
					dA_[i * nx_ + ncols_ + j] = cols->getElements()[rowi];
					posA++;
				}
				else
				{
					dC_[i * nx_ + ncols_ + j] = cols->getElements()[rowi];
					posC++;
				}
			}
		}
	}
	//DSPdebugMessage("posA %d\n", posA);
	for (int i = 0; i < cuts_.sizeCuts(); ++i)
	{
		OsiRowCut * rc = cuts_.rowCutPtr(i);
		if (!rc) continue;
		const CoinPackedVector row = rc->row();
		if (rc->lb() == rc->ub())
		{
			for (int j = 0; j < row.getNumElements(); ++j)
			{
				dA_[posy * nx_ + row.getIndices()[j]] = row.getElements()[j];
				posA++;
			}
			b_[posy] = rc->ub();
			posy++;
		}
		else
		{
			for (int j = 0; j < row.getNumElements(); ++j)
			{
				dC_[posz * nx_ + row.getIndices()[j]] = row.getElements()[j];
				posC++;
			}
			if (rc->lb() > -COIN_DBL_MAX)
			{
				clow_[posz]  = rc->lb();
				iclow_[posz] = 1;
			}
			else
			{
				clow_[posz]  = 0;
				iclow_[posz] = 0;
			}
			if (rc->ub() < COIN_DBL_MAX)
			{
				cupp_[posz]  = rc->ub();
				icupp_[posz] = 1;
			}
			else
			{
				cupp_[posz]  = 0;
				icupp_[posz] = 0;
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
	x_      = new double [nx_];
	y_      = new double [my_];
	lambda_ = new double [mz_];
	pi_     = new double [mz_];
	gamma_  = new double [nx_];
	phi_    = new double [nx_];

	qp_ = new QpGenDense(nx_, my_, mz_);

	/**
	 * Data storage for the regularized master problem:
	 *   Q is the lower triangular matrix of the quadratic term of the objective function.
	 *   A is the equality constraint matrix.
	 *   C is the inequality constraint matrix.
	 */
//	printf("nnzA %d\n", nnzA);
//	for (int i = 0; i < nnzA; ++i)
//		printf(" irowA %d jcolA %d dA %f\n", irowA[i], jcolA[i], dA[i]);

	if (dQ_ == NULL)
	{
		dQ_ = new double [nx_ * nx_];
		CoinZeroN(dQ_, nx_ * nx_);
	}
	prob_ = dynamic_cast<QpGenDense*>(qp_)->makeData(
			c_, dQ_,
			xlow_,  ixlow_,  xupp_,  ixupp_,
			dA_, b_, dC_,
			clow_,  iclow_,  cupp_,  icupp_);
//	prob_->print();

	/** declare variables */
	vars_ = (QpGenVars*)qp_->makeVariables(prob_);

	/** declare residuals */
	resid_ = (QpGenResiduals*)qp_->makeResiduals(prob_);

	/** create solver */
	//solver_ = new GondzioSolver(qp_, prob_);
	solver_ = new MehrotraSolver(qp_, prob_);

	released_ = false;

	END_TRY_CATCH(;)
}

void SolverInterfaceOoqp::gutsOfLoadSparseProblem()
{
	int ndynCols = 0;
	CoinPackedMatrix * cols = NULL; /**< do not free */

	BGN_TRY_CATCH

	/** retrieve dynamic column information */
	if (dynCols_)
	{
		ndynCols = dynCols_->size();
		cols = dynCols_->cols_;
	}

	/** number of variables */
	nx_ = ncols_ + ndynCols;

	/** count number of equality/inequality constraints */
	/** count number of nonzero elements in A and C */
	my_ = 0;
	mz_ = 0;
	nnzA_ = 0;
	nnzC_ = 0;
	for (int i = 0; i < mat_->getNumRows(); ++i)
	{
		if (fabs(rlbd_[i] - rubd_[i]) < 1.0e-10)
		{
			my_++;
			nnzA_ += mat_->getVectorSize(i);
		}
		else
		{
			mz_++;
			nnzC_ += mat_->getVectorSize(i);
		}
	}
	for (int i = 0; i < cuts_.sizeCuts(); ++i)
	{
		OsiRowCut * rc = cuts_.rowCutPtr(i);
		if (!rc) continue;
		if (rc->lb() == rc->ub())
		{
			my_++;
			nnzA_ += rc->row().getNumElements();
		}
		else
		{
			mz_++;
			nnzC_ += rc->row().getNumElements();
		}
	}
	/** for dynamic columns */
	for (int j = 0; j < ndynCols; ++j)
	{
		for (int i = cols->getVectorFirst(j); i < cols->getVectorLast(j); ++i)
		{
			int rowi = cols->getIndices()[i];
			if (fabs(rlbd_[rowi] - rubd_[rowi]) < 1.0e-10)
			{
				//DSPdebugMessage("rowi %d j %d\n", rowi, j);
				nnzA_++;
			}
			else
				nnzC_++;
		}
	}
	//DSPdebugMessage("nnzA %d\n", nnzA);

	/** linear objective */
	FREE_ARRAY_PTR(c_);
	c_ = new double [nx_];
	for (int j = 0; j < ncols_; ++j)
		c_[j] = obj_[j] * sense_;
	/** for dynamic columns */
	for (int j = ncols_; j < nx_; ++j)
		c_[j] = dynCols_->objs_[j-ncols_] * sense_;

	/** column bounds */
	FREE_ARRAY_PTR(xlow_);
	FREE_ARRAY_PTR(ixlow_);
	FREE_ARRAY_PTR(xupp_);
	FREE_ARRAY_PTR(ixupp_);
	xlow_  = new double [nx_];
	ixlow_ = new char [nx_];
	xupp_  = new double [nx_];
	ixupp_ = new char [nx_];
	for (int j = 0; j < ncols_; ++j)
	{
		if (clbd_[j] > -COIN_DBL_MAX)
		{
			xlow_[j] = clbd_[j];
			ixlow_[j] = 1;
		}
		else
		{
			xlow_[j] = 0;
			ixlow_[j] = 0;
		}
		if (cubd_[j] < COIN_DBL_MAX)
		{
			xupp_[j] = cubd_[j];
			ixupp_[j] = 1;
		}
		else
		{
			xupp_[j] = 0;
			ixupp_[j] = 0;
		}
	}

	/** for dynamic columns */
	for (int j = ncols_; j < nx_; ++j)
	{
		if (dynCols_->lbs_[j-ncols_] > -COIN_DBL_MAX)
		{
			xlow_[j] = dynCols_->lbs_[j-ncols_];
			ixlow_[j] = 1;
		}
		else
		{
			xlow_[j] = 0;
			ixlow_[j] = 0;
		}
		if (dynCols_->ubs_[j-ncols_] < COIN_DBL_MAX)
		{
			xupp_[j] = dynCols_->ubs_[j-ncols_];
			ixupp_[j] = 1;
		}
		else
		{
			xupp_[j] = 0;
			ixupp_[j] = 0;
		}
	}

	/** for constraints */
	FREE_ARRAY_PTR(irowA_);
	FREE_ARRAY_PTR(jcolA_);
	FREE_ARRAY_PTR(dA_);
	FREE_ARRAY_PTR(b_);
	FREE_ARRAY_PTR(irowC_);
	FREE_ARRAY_PTR(jcolC_);
	FREE_ARRAY_PTR(dC_);
	FREE_ARRAY_PTR(clow_);
	FREE_ARRAY_PTR(iclow_);
	FREE_ARRAY_PTR(cupp_);
	FREE_ARRAY_PTR(icupp_);
	irowA_ = new int [nnzA_];
	jcolA_ = new int [nnzA_];
	dA_    = new double [nnzA_];
	b_     = new double [my_];
	irowC_ = new int [nnzC_];
	jcolC_ = new int [nnzC_];
	dC_    = new double [nnzC_];
	clow_  = new double [mz_];
	iclow_ = new char [mz_];
	cupp_  = new double [mz_];
	icupp_ = new char [mz_];
	int posA = 0, posC = 0, posy = 0, posz = 0;
	for (int i = 0; i < mat_->getNumRows(); ++i)
	{
		if (fabs(rlbd_[i] - rubd_[i]) < 1.0e-10)
		{
			for (int j = 0; j < mat_->getVectorSize(i); ++j)
			{
				irowA_[posA] = posy;
				jcolA_[posA] = mat_->getIndices()[mat_->getVectorFirst(i) + j];
				dA_[posA] = mat_->getElements()[mat_->getVectorFirst(i) + j];
				posA++;
			}
			b_[posy] = rubd_[i];
			posy++;
		}
		else
		{
			for (int j = 0; j < mat_->getVectorSize(i); ++j)
			{
				irowC_[posC] = posz;
				jcolC_[posC] = mat_->getIndices()[mat_->getVectorFirst(i) + j];
				dC_[posC] = mat_->getElements()[mat_->getVectorFirst(i) + j];
				posC++;
			}
			if (rlbd_[i] > -COIN_DBL_MAX)
			{
				clow_[posz]  = rlbd_[i];
				iclow_[posz] = 1;
			}
			else
			{
				clow_[posz]  = 0;
				iclow_[posz] = 0;
			}
			if (rubd_[i] < COIN_DBL_MAX)
			{
				cupp_[posz]  = rubd_[i];
				icupp_[posz] = 1;
			}
			else
			{
				cupp_[posz]  = 0;
				icupp_[posz] = 0;
			}
			posz++;
		}
		/** for dynamic columns */
		for (int j = 0; j < ndynCols; ++j)
		{
			for (int rowi = cols->getVectorFirst(j); rowi < cols->getVectorLast(j); ++rowi)
			{
				if (i != cols->getIndices()[rowi])
					continue;
				//DSPdebugMessage("irowA %d jcolA %d dA %e\n", i, ncols_ + j, cols->getElements()[rowi]);
				if (fabs(rlbd_[i] - rubd_[i]) < 1.0e-10)
				{
					irowA_[posA] = i;
					jcolA_[posA] = ncols_ + j;
					dA_[posA] = cols->getElements()[rowi];
					posA++;
				}
				else
				{
					irowC_[posC] = i;
					jcolC_[posC] = ncols_ + j;
					dC_[posC] = cols->getElements()[rowi];
					posC++;
				}
			}
		}
	}
	//DSPdebugMessage("posA %d\n", posA);
	for (int i = 0; i < cuts_.sizeCuts(); ++i)
	{
		OsiRowCut * rc = cuts_.rowCutPtr(i);
		if (!rc) continue;
		const CoinPackedVector row = rc->row();
		if (rc->lb() == rc->ub())
		{
			for (int j = 0; j < row.getNumElements(); ++j)
			{
				irowA_[posA] = posy;
				jcolA_[posA] = row.getIndices()[j];
				dA_[posA] = row.getElements()[j];
				posA++;
			}
			b_[posy] = rc->ub();
			posy++;
		}
		else
		{
			for (int j = 0; j < row.getNumElements(); ++j)
			{
				irowC_[posC] = posz;
				jcolC_[posC] = row.getIndices()[j];
				dC_[posC] = row.getElements()[j];
				posC++;
			}
			if (rc->lb() > -COIN_DBL_MAX)
			{
				clow_[posz]  = rc->lb();
				iclow_[posz] = 1;
			}
			else
			{
				clow_[posz]  = 0;
				iclow_[posz] = 0;
			}
			if (rc->ub() < COIN_DBL_MAX)
			{
				cupp_[posz]  = rc->ub();
				icupp_[posz] = 1;
			}
			else
			{
				cupp_[posz]  = 0;
				icupp_[posz] = 0;
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
	x_      = new double [nx_];
	y_      = new double [my_];
	lambda_ = new double [mz_];
	pi_     = new double [mz_];
	gamma_  = new double [nx_];
	phi_    = new double [nx_];

#ifdef DSP_HAS_MA57
	qp_ = new QpGenSparseMa57(nx_, my_, mz_, nnzQ_, nnzA_, nnzC_);
#else
	qp_ = new QpGenSparseMa27(nx_, my_, mz_, nnzQ_, nnzA_, nnzC_);
#endif

	/**
	 * Data storage for the regularized master problem:
	 *   Q is the lower triangular matrix of the quadratic term of the objective function.
	 *   A is the equality constraint matrix.
	 *   C is the inequality constraint matrix.
	 */
//	printf("nnzA %d\n", nnzA);
//	for (int i = 0; i < nnzA; ++i)
//		printf(" irowA %d jcolA %d dA %f\n", irowA[i], jcolA[i], dA[i]);

	prob_ = (QpGenData*)dynamic_cast<QpGenSparseSeq*>(qp_)->copyDataFromSparseTriple(
			c_,     irowQ_, nnzQ_, jcolQ_, dQ_,
			xlow_,  ixlow_,  xupp_,  ixupp_,
			irowA_, nnzA_,   jcolA_, dA_,     b_,
			irowC_, nnzC_,   jcolC_, dC_,
			clow_,  iclow_,  cupp_,  icupp_);
	//prob_->print();

	/** declare variables */
	vars_ = (QpGenVars*)qp_->makeVariables(prob_);

	/** declare residuals */
	resid_ = (QpGenResiduals*)qp_->makeResiduals(prob_);

	/** create solver */
	//solver_ = new GondzioSolver(qp_, prob_);
	solver_ = new MehrotraSolver(qp_, prob_);

	released_ = false;

	END_TRY_CATCH(;)
}

/** solve */
void SolverInterfaceOoqp::solve()
{
	if (released_)
		GUTS_OF_LOAD_PROBLEM();

	//prob_->print();
	if (par_->getIntParam("LOG_LEVEL") >= 10)
		solver_->monitorSelf();
	int status = solver_->solve(prob_, vars_, resid_);
	if (par_->getIntParam("LOG_LEVEL") >= 10)
		printf("OOQP status: %d\n", status);
	switch(status)
	{
	case 0:
	{
		status_ = DSP_STAT_OPTIMAL;
//		double phi = fabs(resid_->residualNorm() - resid_->dualityGap()) / prob_->datanorm();
//		if (phi < 1.0e-6)
//			status_ = STO_STAT_OPTIMAL;
//		else
//			status_ = STO_STAT_UNKNOWN;
		objval_ = prob_->objectiveValue(vars_);
		vars_->x->copyIntoArray(x_);
		vars_->y->copyIntoArray(y_);
		vars_->lambda->copyIntoArray(lambda_);
		nLambdas_ =  vars_->lambda->length();
		vars_->pi->copyIntoArray(pi_);
		nPis_ = vars_->pi->length();
		vars_->gamma->copyIntoArray(gamma_);
		vars_->phi->copyIntoArray(phi_);
		nIters_ = solver_->iter;
		break;
	}
	case 1:
		status_ = DSP_STAT_STOPPED_UNKNOWN;
		break;
	case 2:
		status_ = DSP_STAT_STOPPED_ITER;
		break;
	case 3:
		status_ = DSP_STAT_PRIM_INFEASIBLE;
		break;
	case 4:
	default:
		status_ = DSP_STAT_UNKNOWN;
		break;
	}
}

/** add row */
void SolverInterfaceOoqp::addRow(int size, const int * indices, const double * vals, double lb, double ub)
{
	releaseOOQP();
	/** TODO: wrong */
	OsiRowCut rc;
	rc.setRow(size, indices, vals);
	rc.setLb(lb);
	rc.setUb(ub);
	cuts_.insert(rc);
}

/** add column */
void SolverInterfaceOoqp::addCol(int size, const int * indices, const double * vals, double lb, double ub, double obj)
{
	assert(cuts_.sizeCuts() == 0);

	if (dynCols_ == NULL)
		dynCols_ = new dynColumns(nrows_, ncols_);

	dynCols_->addCol(size, indices, vals, lb, ub, obj);
}

/** clear columns */
void SolverInterfaceOoqp::clearCols()
{
	if (dynCols_ != NULL) dynCols_->clear();
}

/** add cuts */
int SolverInterfaceOoqp::addCuts(OsiCuts cuts, double effectiveness)
{
	releaseOOQP();
	cuts_.insert(cuts);
	return cuts.sizeCuts();
}

/** set objective coefficients */
void SolverInterfaceOoqp::setObjCoef(double * obj)
{
	releaseOOQP();
	CoinCopyN(obj, ncols_, obj_);
}

/** set lower triangular Hessian matrix */
void SolverInterfaceOoqp::setHessian(int nnzQ, int * irowQ, int * jcolQ, double * dQ)
{
	SET_HESSIAN(nnzQ, irowQ, jcolQ, dQ);
}

/** set lower triangular Hessian matrix */
void SolverInterfaceOoqp::setSparseHessian(int nnzQ, int * irowQ, int * jcolQ, double * dQ)
{
	releaseOOQP();
	FREE_ARRAY_PTR(irowQ_);
	FREE_ARRAY_PTR(jcolQ_);
	FREE_ARRAY_PTR(dQ_);
	nnzQ_ = nnzQ;
	if (nnzQ_ == 0)
	{
		irowQ_ = NULL;
		jcolQ_ = NULL;
		dQ_ = NULL;
	}
	else
	{
		irowQ_ = new int [nnzQ_];
		jcolQ_ = new int [nnzQ_];
		dQ_ = new double [nnzQ_];
		CoinCopyN(irowQ, nnzQ_, irowQ_);
		CoinCopyN(jcolQ, nnzQ_, jcolQ_);
		for (int i = 0; i <nnzQ_; ++i)
			dQ_[i] = dQ[i] * sense_;
	}
}

/** set lower triangular Hessian matrix */
void SolverInterfaceOoqp::setDenseHessian(int nnzQ, int * irowQ, int * jcolQ, double * dQ)
{
	releaseOOQP();
	FREE_ARRAY_PTR(dQ_);
	nnzQ_ = nx_ * nx_;
	dQ_ = new double [nnzQ_];
	CoinZeroN(dQ_, nnzQ_);
	for (int i = 0; i < nnzQ; ++i)
		dQ_[irowQ[i]*nx_+jcolQ[i]] = dQ[i] * sense_;
}

/** set column bounds */
void SolverInterfaceOoqp::setColLower(int index, double lb)
{
	releaseOOQP();
	clbd_[index] = lb;
}

/** set column bounds */
void SolverInterfaceOoqp::setColUpper(int index, double ub)
{
	releaseOOQP();
	cubd_[index] = ub;
}

/** set column bounds */
void SolverInterfaceOoqp::setColBounds(int index, double lb, double ub)
{
	releaseOOQP();
	clbd_[index] = lb;
	cubd_[index] = ub;
}

/** set row bounds */
void SolverInterfaceOoqp::setRowLower(int index, double lb)
{
	releaseOOQP();
	rlbd_[index] = lb;
}

/** set row bounds */
void SolverInterfaceOoqp::setRowUpper(int index, double ub)
{
	releaseOOQP();
	rubd_[index] = ub;
}

/** set row bounds */
void SolverInterfaceOoqp::setRowBounds(int index, double lb, double ub)
{
	releaseOOQP();
	rlbd_[index] = lb;
	rubd_[index] = ub;
}
