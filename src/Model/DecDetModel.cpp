/*
 * DecDetModel.cpp
 *
 *  Created on: Jul 1, 2015
 *      Author: ctjandra
 */

#include "DecDetModel.h"

/** 
 * Decomposable version of a deterministic model.
 */

/** copy constructor */
DecDetModel::DecDetModel(const DecDetModel & rhs)
{
	decData_ = new DecData(*rhs.decData_);
	det_ = new DetModel(*rhs.det_);
}

DecDetModel::~DecDetModel()
{
	FREE_PTR(decData_);
	FREE_PTR(det_);
}

STO_RTN_CODE DecDetModel::decompose(
		int size,                /**< [in] size of subproblem subset */
		int * subprobs,          /**< [in] subset of subproblems */
		int naux,                /**< [in] number of auxiliary columns */
		double * clbd_aux,       /**< [in] lower bounds for auxiliary columns */
		double * cubd_aux,       /**< [in] upper bounds for auxiliary columns */
		double * obj_aux,        /**< [in] objective coefficients for auxiliary columns */
		CoinPackedMatrix *& mat, /**< [out] constraint matrix */
		double *& clbd,          /**< [out] column lower bounds */
		double *& cubd,          /**< [out] column upper bounds */
		char   *& ctype,         /**< [out] column types */
		double *& obj,           /**< [out] objective coefficients */
		double *& rlbd,          /**< [out] row lower bounds */
		double *& rubd           /**< [out] row upper bounds */)
{
	BGN_TRY_CATCH

	int nrows = det_->getNumRows();
	int ncols = det_->getNumCols();

	/** allocate memory */
	clbd = new double [ncols + naux];
	cubd = new double [ncols + naux];
	ctype = new char [ncols + naux];
	obj = new double [ncols + naux];
	rlbd = new double [nrows];
	rubd = new double [nrows];

	/** copy data from DetModel */
	mat = new CoinPackedMatrix(*(det_->getConstraintMatrix()));
	CoinCopyN(det_->getColLower(), ncols, clbd);
	CoinCopyN(det_->getColUpper(), ncols, cubd);
	CoinCopyN(det_->getCtype(),    ncols, ctype);
	CoinCopyN(det_->getObj(),      ncols, obj);
	CoinCopyN(det_->getRowLower(), nrows, rlbd);
	CoinCopyN(det_->getRowUpper(), nrows, rubd);

	decData_->decompose(size, subprobs, mat, clbd, cubd, ctype, obj, rlbd, rubd);

	/** auxiliary columns */
	CoinCopyN(clbd_aux, naux, clbd + ncols);
	CoinCopyN(cubd_aux, naux, cubd + ncols);
	CoinCopyN(obj_aux, naux, obj + ncols);
	CoinFillN(ctype + ncols, naux, 'C');

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
}

STO_RTN_CODE DecDetModel::getFullModel(
	CoinPackedMatrix *& mat, /**< [out] constraint matrix */
	double *& clbd,          /**< [out] column lower bounds */
	double *& cubd,          /**< [out] column upper bounds */
	char   *& ctype,         /**< [out] column types */
	double *& obj,           /**< [out] objective coefficients */
	double *& rlbd,          /**< [out] row lower bounds */
	double *& rubd           /**< [out] row upper bounds */)
{
	BGN_TRY_CATCH

	int nrows = det_->getNumRows();
	int ncols = det_->getNumCols();
	int nrows_coupling = decData_->getNumCouplingRows();

	/** allocate memory */
	clbd = new double [ncols];
	cubd = new double [ncols];
	ctype = new char [ncols];
	obj = new double [ncols];
	rlbd = new double [nrows + nrows_coupling];
	rubd = new double [nrows + nrows_coupling];

	/** copy data from DetModel */
	assert(det_->getConstraintMatrix()->getNumCols() == decData_->getConstraintMatrix()->getNumCols());
	mat = new CoinPackedMatrix(*(det_->getConstraintMatrix()));
	mat->bottomAppendPackedMatrix(*(decData_->getConstraintMatrix()));
	CoinCopyN(det_->getColLower(), ncols, clbd);
	CoinCopyN(det_->getColUpper(), ncols, cubd);
	CoinCopyN(det_->getCtype(),    ncols, ctype);
	CoinCopyN(det_->getObj(),      ncols, obj);
	CoinCopyN(det_->getRowLower(), nrows, rlbd);
	CoinCopyN(det_->getRowUpper(), nrows, rubd);

	/** TODO: Only equality to zero is supported for now */
	CoinZeroN(rlbd + nrows, nrows_coupling);
	CoinZeroN(rubd + nrows, nrows_coupling);

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
}

/** get objective coefficients */
void DecDetModel::getObjCoef(double * obj)
{
	CoinCopyN(det_->getObj(), det_->getNumCols(), obj);
}
