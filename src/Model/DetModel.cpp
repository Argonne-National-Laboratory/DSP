/*
 * DetModel.cpp
 *
 *  Created on: July 1, 2015
 *      Author: ctjandra
 */

#include "DetModel.h"

DetModel::DetModel() :
		mat_(NULL),
		clbd_(NULL),
		cubd_(NULL),
		ctype_(NULL),
		obj_(NULL),
		rlbd_(NULL),
		rubd_(NULL),
		nints_(0)
{
	/** nothing to do */
}

DetModel::DetModel(
	const CoinPackedMatrix * mat,   /**< constraint matrix */
	const double * clbd,            /**< column lower bounds */
	const double * cubd,            /**< column upper bounds */
	const char   * ctype,           /**< column types */
	const double * obj,             /**< objective coefficients */
	const double * rlbd,            /**< row lower bounds */
	const double * rubd             /**< row upper bounds */)
{
	int ncols = mat->getNumCols();
	int nrows = mat->getNumRows();

	clbd_ = new double [ncols];
	cubd_ = new double [ncols];
	ctype_ = new char [ncols];
	obj_ = new double [ncols];
	rlbd_ = new double [nrows];
	rubd_ = new double [nrows];

	mat_ = new CoinPackedMatrix(*mat);
	CoinCopyN(clbd,  ncols, clbd_);
	CoinCopyN(cubd,  ncols, cubd_);
	CoinCopyN(ctype, ncols, ctype_);
	CoinCopyN(obj,   ncols, obj_);
	CoinCopyN(rlbd,  nrows, rlbd_);
	CoinCopyN(rubd,  nrows, rubd_);

	nints_ = 0;
	for (int j = 0; j < ncols; j++)
	{
		if (ctype_[j] != 'C')
			nints_++;
	}
}

DetModel::DetModel(
	const CoinBigIndex * start,     /**< start index for each row */
	const int    * index,           /**< column indices */
	const double * value,           /**< constraint elements */
	const int      numels,          /**< number of elements in index and value */
	const int      ncols,           /**< number of columns */
	const int      nrows,           /**< number of rows */
	const double * clbd,            /**< column lower bounds */
	const double * cubd,            /**< column upper bounds */
	const char   * ctype,           /**< column types */
	const double * obj,             /**< objective coefficients */
	const double * rlbd,            /**< row lower bounds */
	const double * rubd             /**< row upper bounds */)
{
	int * len = new int [nrows];
	for (int i = 0; i < nrows; i++)
		len[i] = start[i+1] - start[i];

	mat_ = new CoinPackedMatrix(false, ncols, nrows, numels, value, index, start, len);

	clbd_ = new double [ncols];
	cubd_ = new double [ncols];
	ctype_ = new char [ncols];
	obj_ = new double [ncols];
	rlbd_ = new double [nrows];
	rubd_ = new double [nrows];

	CoinCopyN(clbd,  ncols, clbd_);
	CoinCopyN(cubd,  ncols, cubd_);
	CoinCopyN(ctype, ncols, ctype_);
	CoinCopyN(obj,   ncols, obj_);
	CoinCopyN(rlbd,  nrows, rlbd_);
	CoinCopyN(rubd,  nrows, rubd_);

	nints_ = 0;
	for (int j = 0; j < ncols; j++)
	{
		if (ctype_[j] != 'C')
			nints_++;
	}
}

/** copy constructor */
DetModel::DetModel(const DetModel & rhs) :
		nints_(rhs.nints_)
{
	int ncols = rhs.mat_->getNumCols();
	int nrows = rhs.mat_->getNumRows();

	clbd_ = new double [ncols];
	cubd_ = new double [ncols];
	ctype_ = new char [ncols];
	obj_ = new double [ncols];
	rlbd_ = new double [nrows];
	rubd_ = new double [nrows];

	mat_ = new CoinPackedMatrix(*rhs.mat_);
	CoinCopyN(rhs.clbd_,  ncols, clbd_);
	CoinCopyN(rhs.cubd_,  ncols, cubd_);
	CoinCopyN(rhs.ctype_, ncols, ctype_);
	CoinCopyN(rhs.obj_,   ncols, obj_);
	CoinCopyN(rhs.rlbd_,  nrows, rlbd_);
	CoinCopyN(rhs.rubd_,  nrows, rubd_);
}

DetModel::~DetModel()
{
	FREE_PTR(mat_);
	FREE_ARRAY_PTR(clbd_);
	FREE_ARRAY_PTR(cubd_);
	FREE_ARRAY_PTR(ctype_);
	FREE_ARRAY_PTR(obj_);
	FREE_ARRAY_PTR(rlbd_);
	FREE_ARRAY_PTR(rubd_);
	nints_ = 0;
}

void DetModel::__printData()
{
	int ncols = mat_->getNumCols();
	int nrows = mat_->getNumRows();

	printf("\n### BEGINNING of printing DetModel data ###\n\n");

	printf("nrows %d\n", nrows);
	printf("ncols %d\n", ncols);
	PRINT_ARRAY_MSG(ncols, clbd_, "clbd_")
	PRINT_ARRAY_MSG(ncols, cubd_, "cubd_")
	PRINT_ARRAY_MSG(ncols, ctype_, "ctype_")
	PRINT_ARRAY_MSG(ncols, obj_, "obj_")
	PRINT_ARRAY_MSG(nrows, rlbd_, "rlbd_")
	PRINT_ARRAY_MSG(nrows, rubd_, "rubd_")

	printf("=== BEGINNING of CoinPackedMatrix mat_ ===\n");
	printf("isColOrdered %d\n", mat_->isColOrdered());
	PRINT_ARRAY_MSG(mat_->getMajorDim(), mat_->getVectorStarts(), "VectorStarts")
	PRINT_SPARSE_ARRAY_MSG(
			mat_->getNumElements(),
			mat_->getIndices(),
			mat_->getElements(),
			"Elements")
	printf("=== END of CoinPackedMatrix mat_ ===\n");

	printf("\n### END of printing StoModel data ###\n");
}
