/*
 * DecTssModel.cpp
 *
 *  Created on: Jun 17, 2015
 *      Author: ctjandra
 */

#include "DecTssModel.h"

DecTssModel::DecTssModel() :
TssModel()
{
	/** nothing to do */
}

/** copy constructor */
DecTssModel::DecTssModel(const DecTssModel & rhs) :
TssModel(rhs)
{
	/** nothing to do */
}

/** copy constructor */
DecTssModel::DecTssModel(const TssModel & rhs) :
TssModel(rhs)
{
	/** nothing to do */
}
DecTssModel::~DecTssModel()
{
	/** nothing to do */
}

DSP_RTN_CODE DecTssModel::decomposeCoupling(
	int size,                    /**< [in] size of subproblem subset */
	int * subprobs,              /**< [in] subset of subproblems */
	CoinPackedMatrix *& cpl_mat, /**< [out] coupling constraint matrix */
	int *& cpl_cols,             /**< [out] columns of cpl_mat involved in coupling rows */
	int & cpl_ncols              /**< [out] size of cpl_cols */)
{
	int ncols_first = getNumCols(0);
	int ncols = getNumCols(0) + size * getNumCols(1);

	/** cpl_mat is the identity matrix of size ncols: nonanticipativity is in form x = x_s,
	  * but with x implicitly eliminated, we only have x_s in each row. */
	double * elem = new double [ncols_first];
	int * ind = new int [ncols_first];
	int * start = new int [ncols_first+1];
	int * len = new int [ncols_first];
	for (int i = 0; i < ncols_first; i++)
	{
		ind[i] = i;
		start[i] = i;
		len[i] = 1;
		elem[i] = 1.0;
	}
	start[ncols_first] = ncols_first;
	cpl_mat = new CoinPackedMatrix(false, ncols, ncols_first, ncols_first, elem, ind, start, len);

	/** columns involved are the first ncols_first columns */
	cpl_ncols = ncols_first;
	cpl_cols = new int [cpl_ncols];
	for (int i = 0; i < cpl_ncols; i++)
		cpl_cols[i] = i;

	int k = 0;
	for (int s = 0; s < size; s++)
		k += getNumSubproblemCouplingCols(subprobs[s]);
	assert(k == cpl_ncols);

	return DSP_RTN_OK;
}

const int * DecTssModel::getSubproblemCouplingRowIndices(int s)
{
	int * rowIndices = new int [getNumCols(0)];
	for (int i = 0; i < getNumCols(0); i++)
		rowIndices[i] = s * getNumCols(0) + i;
	return rowIndices;
}

double DecTssModel::evalLhsCouplingRow(int row, double ** solutions)
{
	double val = 0;
	for (int s = 0; s < getNumSubproblems(); s++)
		val += evalLhsCouplingRowSubprob(row, s, solutions[s]);
	return val;
}

double DecTssModel::evalLhsCouplingRowSubprob(int row, int s, double * subprobSolution)
{
	/** recall solution is in the space of the coupling vars. of subproblem,
	 *  while row is from the master problem coupling matrix */
	if (row < s * getNumCols(0) || row >= (s + 1) * getNumCols(0))
		return 0;
	double sol = subprobSolution[row - s * getNumCols(0)];
	if (fabs(sol) <= 1E-10)
		return 0;
	return sol;
}

double DecTssModel::evalLhsCouplingRowAlternative(int row, double ** solutions)
{
	double val = 0;
//	for (int s = 0; s < getNumSubproblems(); s++)
//		if (row + 1 < getNumCols(0))
//			val += solutions[s][row] - solutions[s][row+1];
//		else
//			val += solutions[s][row] - solutions[s][0];
	int s = row / getNumSubproblems();
	int i = row % getNumSubproblems();
	val = solutions[s][i];
	if (s + 1 < getNumSubproblems())
		val -= solutions[s+1][i];
	else
		val -= solutions[0][i];
	return val;
}

void DecTssModel::convertLagrangianFromAlternative(double * multipliers, double *& newMultipliers)
{
	int ncols_first = getNumCols(0);
	int nsubprobs = getNumSubproblems();
	for (int s = 0; s < nsubprobs; ++s)
		for (int j = 0; j < ncols_first; ++j)
		{
			if (s == 0)
				newMultipliers[s * ncols_first + j] = multipliers[s * ncols_first + j] - multipliers[(nsubprobs - 1) * ncols_first + j];
			else
				newMultipliers[s * ncols_first + j] = multipliers[s * ncols_first + j] - multipliers[(s - 1) * ncols_first + j];
		}
}

DSP_RTN_CODE DecTssModel::getFullModel(
	CoinPackedMatrix *& mat, /**< [out] constraint matrix */
	double *& clbd,          /**< [out] column lower bounds */
	double *& cubd,          /**< [out] column upper bounds */
	char   *& ctype,         /**< [out] column types */
	double *& obj,           /**< [out] objective coefficients */
	double *& rlbd,          /**< [out] row lower bounds */
	double *& rubd           /**< [out] row upper bounds */)
{
	int nsubprobs = getNumSubproblems();
	int * subprobs = new int [nsubprobs];
	CoinIotaN(subprobs, nsubprobs, 0);

	DSP_RTN_CODE rtn = decompose(nsubprobs, subprobs, 0, NULL, NULL, NULL,
			mat, clbd, cubd, ctype, obj, rlbd, rubd);

	/** free array */
	FREE_ARRAY_PTR(subprobs);

	return rtn;
}

/** get objective coefficients */
void DecTssModel::getObjCoef(double * obj)
{
	CoinCopyN(obj_core_[0], ncols_[0], obj);
	for (int s = 0; s < nscen_; ++s)
	{
		CoinCopyN(obj_core_[1], ncols_[1], obj + ncols_[0] + s * ncols_[1]);
		for (int j = 0; j < obj_scen_[s]->getNumElements(); ++j)
		{
			obj[s * ncols_[1] + obj_scen_[s]->getIndices()[j]] = obj_scen_[s]->getElements()[j];
		}
		for (int j = 0; j < ncols_[1]; ++j)
			obj[ncols_[0] + s * ncols_[1] + j] *= prob_[s];
	}
}
