/*
 * DecData.cpp
 *
 *  Created on: Jun 17, 2015
 *      Author: ctjandra
 */

/** 
 * Data structure to explicitly store coupling constraints and column partition, meant to be attached
 * to a decomposition model.
 */

#include "DecData.h"

/**
 * Used to check if a coefficient is nonzero: we should not have to worry about floating point precision
 * here because the matrices should come directly from the model, but this may be changed if we do.
 */
#define NONZERO(x) ((x) != 0)

/**
 * Solution variables with at most this value (in absolute terms) are considered to be zero for purposes
 * of evaluating it on coupling constraints.
 */
#define EPSILON_ROW_EVAL 1E-10

/** copy constructor */
DecData::DecData(const DecData & rhs) :
		nsubprobs_(rhs.nsubprobs_),
		nrows_coupling_(rhs.nrows_coupling_),
		ncols_coupling_(rhs.ncols_coupling_),
		ncols_(rhs.ncols_)
{
	cols_subprob_ = new int [ncols_];
	ncols_subprob_ = new int [nsubprobs_];
	ncols_subprob_coupling_ = new int [nsubprobs_];
	nrows_subprob_coupling_ = new int [nsubprobs_];
	cols_incoupling_ = new bool [ncols_];
	rows_subprob_ = new int * [nsubprobs_];
	orig_coupling_col_ = new int * [nsubprobs_];

	CoinCopyN(rhs.cols_subprob_,           ncols_,     cols_subprob_);
	CoinCopyN(rhs.ncols_subprob_,          nsubprobs_, ncols_subprob_);
	CoinCopyN(rhs.ncols_subprob_coupling_, nsubprobs_, ncols_subprob_coupling_);
	CoinCopyN(rhs.nrows_subprob_coupling_, nsubprobs_, nrows_subprob_coupling_);
	CoinCopyN(rhs.cols_incoupling_,        ncols_,     cols_incoupling_);

	for (int s = 0; s < nsubprobs_; s++)
	{
		orig_coupling_col_[s] = new int [ncols_];
		rows_subprob_[s] = new int [nrows_coupling_];
		CoinCopyN(rhs.orig_coupling_col_[s], ncols_, orig_coupling_col_[s]);
		CoinCopyN(rhs.rows_subprob_[s], nrows_coupling_, rows_subprob_[s]);
	}

	rows_coupling_ = new CoinPackedMatrix(*rhs.rows_coupling_);
}

DecData::~DecData()
{
	FREE_ARRAY_PTR(nrows_subprob_coupling_);
	FREE_ARRAY_PTR(ncols_subprob_coupling_);
	FREE_ARRAY_PTR(ncols_subprob_);
	FREE_ARRAY_PTR(cols_subprob_);
	FREE_ARRAY_PTR(cols_incoupling_);
	FREE_2D_ARRAY_PTR(nsubprobs_, rows_subprob_);
	FREE_2D_ARRAY_PTR(nsubprobs_, orig_coupling_col_);
	FREE_PTR(rows_coupling_);
}

DecData::DecData(
		int      nsubprobs,      /**< number of subproblems */
		int      ncols,          /**< number of columns for the entire problem */
		int      ncoupling,      /**< number of coupling constraints */
		int *    varPartition,   /**< partition of columns into subproblems */
		int *    couplingStarts, /**< indices in cols at which each coupling constraint starts */
		int *    couplingCols,   /**< variables of each coupling constraint left-hand side */
		double * couplingCoeffs  /**< coefficients of each coupling constraint left-hand side */)
{
	ncols_ = ncols;
	nsubprobs_ = nsubprobs;
	nrows_coupling_ = ncoupling;
	cols_subprob_ = new int [ncols_];
	for (int c = 0; c < ncols_; c++)
		cols_subprob_[c] = varPartition[c];

	ncols_subprob_ = new int [nsubprobs_];
	ncols_subprob_coupling_ = new int [nsubprobs_];
	nrows_subprob_coupling_ = new int [nsubprobs_];
	for (int s = 0; s < nsubprobs_; s++)
	{
		ncols_subprob_[s] = 0;
		ncols_subprob_coupling_[s] = 0;
		nrows_subprob_coupling_[s] = 0;
	}

	/** construct coupling constraints */
	rows_coupling_ = new CoinPackedMatrix(false, 0, 0);
	rows_coupling_->setDimensions(0, ncols_);
	cols_incoupling_ = new bool [ncols_];
	for (int c = 0; c < ncols_; c++)
		cols_incoupling_[c] = false;
	for (int i = 0; i < nrows_coupling_; i++)
	{
		int size = couplingStarts[i+1] - couplingStarts[i];
		int * inds = new int [size];
		double * elems = new double [size];
		int k = 0;
		for (int j = couplingStarts[i]; j < couplingStarts[i+1]; j++)
		{
			assert(couplingCols[j] >= 0 && couplingCols[j] < ncols);
			cols_incoupling_[couplingCols[j]] = true;
			inds[k] = couplingCols[j];
			elems[k] = couplingCoeffs[j];
			k++;
		}
		rows_coupling_->appendRow(size, inds, elems);
	}

	/** identify and count columns corresponding to subproblems */
	ncols_coupling_ = 0;
	orig_coupling_col_ = new int * [nsubprobs_];
	for (int s = 0; s < nsubprobs_; s++)
		orig_coupling_col_[s] = new int [ncols_];
	for (int i = 0; i < ncols_; i++)
	{
		/** subproblem corresponding to column i */
		int s = cols_subprob_[i];
		assert(s >= 0 && s < nsubprobs);
		if (cols_incoupling_[i])
		{
			orig_coupling_col_[s][ncols_subprob_coupling_[s]] = i;
			ncols_coupling_++;
			ncols_subprob_coupling_[s]++;
		}
		ncols_subprob_[s]++;
	}

	/** identify and count coupling rows corresponding to subproblems */
	rows_subprob_ = new int * [nsubprobs_];
	for (int s = 0; s < nsubprobs_; s++)
		rows_subprob_[s] = new int [nrows_coupling_];
	for (int i = 0; i < nrows_coupling_; i++)
	{
		const int * inds = rows_coupling_->getVector(i).getIndices();
		const double * elems = rows_coupling_->getVector(i).getElements();
		for (int j = 0; j < rows_coupling_->getVector(i).getNumElements(); j++)
		{
			int s = cols_subprob_[inds[j]];
			if (NONZERO(elems[j]))
			{
				rows_subprob_[s][nrows_subprob_coupling_[s]] = i;
				nrows_subprob_coupling_[s]++;
			}
		}
	}
}

double DecData::evalLhsRow(int row, double ** solutions)
{
	assert(row >= 0 && row < nrows_coupling_);
	assert(!rows_coupling_->isColOrdered());

	double val = 0;
	for (int s = 0; s < nsubprobs_; s++)
	{
		val += evalLhsRowSubprob(row, s, solutions[s]);
	}
	return val;
}

double DecData::evalLhsRowSubprob(int row, int subprob, double * subprobSolution)
{
	assert(row >= 0 && row < nrows_coupling_);
	assert(!rows_coupling_->isColOrdered());

	double val = 0;
	for (int i = 0; i < ncols_subprob_coupling_[subprob]; i++)
	{
		if (fabs(subprobSolution[i]) > EPSILON_ROW_EVAL)
		{
			/** subprobSolution is in the space of the subproblem:
			 *  need to recover column from master problem space */
			int col = orig_coupling_col_[subprob][i];
			val += rows_coupling_->getVector(row)[col] * subprobSolution[i];
		}
	}
	return val;
}

STO_RTN_CODE DecData::decompose(
	int size,                /**< [in] size of subproblem subset */
	int * subprobs,          /**< [in] subset of subproblems */
	CoinPackedMatrix *& mat, /**< [in/out] constraint matrix */
	double *& clbd,          /**< [in/out] column lower bounds */
	double *& cubd,          /**< [in/out] column upper bounds */
	char   *& ctype,         /**< [in/out] column types */
	double *& obj,           /**< [in/out] objective coefficients */
	double *& rlbd,          /**< [in/out] row lower bounds */
	double *& rubd           /**< [in/out] row upper bounds */)
{
	if (ncols_ != mat->getNumCols())
	{
		printf("Error: Matrix supplied for decomposition does not match decomposition specifications\n");
		return STO_RTN_ERR;
	}

	if (!checkPartitionConsistency(mat))
	{
		printf("Error: Partition is inconsistent: there cannot be two nonzero coefficients of different subproblems in the same (non-coupling) row.\n");
		return STO_RTN_ERR;
	}

	int orig_nrows = mat->getNumRows();

	/** keep only rows and columns that are relevant to the subproblems */
	vector<int> rowsDiscarded;
	vector<int> colsDiscarded;
	reduceMatrixToSubproblems(size, subprobs, mat, rowsDiscarded, colsDiscarded);
	assert(mat->getNumCols() + colsDiscarded.size() == ncols_);
	assert(mat->getNumRows() + rowsDiscarded.size() == orig_nrows);

	/** update bounds and coefficients */
	for (int j = 0, i = 0; j + i < ncols_; j++)
	{
		while (i < colsDiscarded.size() && colsDiscarded[i] == j + i)
			i++;
		clbd[j] = clbd[j+i];
		cubd[j] = cubd[j+i];
		ctype[j] = ctype[j+i];
		obj[j] = obj[j+i];
	}

	for (int j = 0, i = 0; j + i < orig_nrows; j++)
	{
		while (i < rowsDiscarded.size() && rowsDiscarded[i] == j + i)
			i++;
		rlbd[j] = rlbd[j+i];
		rubd[j] = rubd[j+i];
	}

	return STO_RTN_OK;
}

STO_RTN_CODE DecData::decomposeCoupling(
		int size,                    /**< [in] size of subproblem subset */
		int * subprobs,              /**< [in] subset of subproblems */
		CoinPackedMatrix *& cpl_mat, /**< [out] coupling constraint matrix */
		int *& cpl_cols,             /**< [out] columns of cpl_mat involved in coupling rows */
		int & cpl_ncols              /**< [out] size of cpl_cols */)
{
	/** allocate memory */
	cpl_mat = new CoinPackedMatrix(*rows_coupling_);
	cpl_cols = new int [ncols_];

	/** decompose coupling constraints */
	vector<int> rowsDiscarded;
	vector<int> colsDiscarded;
	reduceMatrixToSubproblems(size, subprobs, cpl_mat, rowsDiscarded, colsDiscarded);

	/** get coupling column indices for subproblem */
	int * colMap = new int [ncols_]; /**< auxiliary map from original space to subproblem space
									   *  (i.e. the i-th column of rows_coupling_ is the colMap[i]-th column of new cpl_mat) */
	int j = 0;
	for (int c = 0; c < ncols_; c++)
	{
		if (j < colsDiscarded.size() && colsDiscarded[j] == c)
		{
			colMap[c] = -1;
			j++;
		}
		else
			colMap[c] = c - j;
	}

	cpl_ncols = 0;
	for (int c = 0; c < ncols_; c++)
	{
		if (cols_incoupling_[c] && columnInSubprobs(c, size, subprobs))
		{
			assert(colMap[c] >= 0);
			cpl_cols[cpl_ncols++] = colMap[c];
		}
	}

	int k = 0;
	for (int s = 0; s < size; s++)
		k += ncols_subprob_coupling_[subprobs[s]];
	assert(k == cpl_ncols);

	return STO_RTN_OK;
}

bool DecData::columnInSubprobs(
	int col,                 /**< column to check if in subproblem subset */
	int size,                /**< size of subproblem subset */
	int * subprobs           /**< subset of subproblems */)
{
	for (int s = 0; s < size; s++)
		if (cols_subprob_[col] == subprobs[s])
			return true;
	return false;
}

void DecData::reduceMatrixToSubproblems(
	int size,                    /**< [in] size of subproblem subset */
	int * subprobs,              /**< [in] subset of subproblems */
	CoinPackedMatrix *& mat,     /**< [in/out] constraint matrix */
	vector<int> & rowsDiscarded, /**< [out] row indices from the matrix that were discarded */
	vector<int> & colsDiscarded  /**< [out] column indices from the matrix that were discarded */)
{
	/** discard columns that do not correspond to given subproblems */
	colsDiscarded.clear();
	for (int c = 0; c < ncols_; c++)
		if (!columnInSubprobs(c, size, subprobs))
			colsDiscarded.push_back(c);

	mat->deleteCols(colsDiscarded.size(), &colsDiscarded[0]);

	/** discard rows of zeroes (which correspond to other subproblems) */
	rowsDiscarded.clear();
	if (mat->isColOrdered())
	{
		/** column ordered */
		bool * zeroRows = new bool [mat->getNumRows()]; /**< auxiliary vector to flag zero rows */
		for (int i = 0; i < mat->getNumRows(); i++)
			zeroRows[i] = true;

		/** unmark any index (row) with nonzero element */
		for (int i = 0; i < mat->getNumElements(); i++)
			if (NONZERO(mat->getElements()[i]))
				zeroRows[mat->getIndices()[i]] = false;

		for (int i = 0; i < mat->getNumRows(); i++)
			if (zeroRows[i])
				rowsDiscarded.push_back(i);
	}
	else
	{
		/** row ordered */
		for (int i = 0; i < mat->getNumRows(); i++)
		{
			const double * rowElems = mat->getVector(i).getElements();
			bool found = false;
			for (int j = 0; j < mat->getVectorSize(i); j++)
			{
				if (NONZERO(rowElems[j]))
				{
					found = true;
					break;
				}
			}

			if (!found)
				rowsDiscarded.push_back(i);
		}
	}

	mat->deleteRows(rowsDiscarded.size(), &rowsDiscarded[0]);
}

/**
 * Return false if there are two nonzero coefficients of different classes in the same row
 * of a given matrix.
 */
bool DecData::checkPartitionConsistency(const CoinPackedMatrix * mat)
{
	assert(mat->getNumCols() == ncols_); /** matrix should be compatible with partition */
	if (mat->isColOrdered()) /** matrix should be row ordered */
	{
		/* TODO Implement consistency check for column ordered matrices if necessary */
		printf("Warning: Matrices should be row ordered, consistency not checked\n");
		return true; /** refuse to check consistency */
	}

	for (int i = 0; i < mat->getNumRows(); i++)
	{
		int size = mat->getVector(i).getNumElements();
		const int * inds = mat->getVector(i).getIndices();
		const double * elems = mat->getVector(i).getElements();
		int partitionClass = -1;
		for (int j = 0; j < size; j++)
		{
			if (NONZERO(elems[j]))
			{
				if (partitionClass == -1)
					partitionClass = cols_subprob_[inds[j]];
				else if (partitionClass != cols_subprob_[inds[j]])
					return false; /** inconsistent */
			}
		}
	}

	return true;
}
