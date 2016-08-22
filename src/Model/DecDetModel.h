/*
 * DecDetModel.h
 *
 *  Created on: Jul 1, 2015
 *      Author: ctjandra
 */

#ifndef DECDETMODEL_H_
#define DECDETMODEL_H_

#include "TssModel.h"
#include "DecModel.h"
#include "DecData.h"

/** 
 * Decomposable version of a deterministic model.
 */
class DecDetModel: public DecModel {

public:

	/** default constructor */
	DecDetModel() :
		det_(NULL),
		decData_(NULL)
	{
		/** nothing to do */
	}

	DecDetModel(DetModel * det, DecData * decData) :
		det_(det)
	{
		decData_ = NULL;
		if (decData)
			decData_ = new DecData(*decData);
	}

	/** copy constructor */
	DecDetModel(const DecDetModel & rhs);

	/** default destructor */
	virtual ~DecDetModel();

public:

	DecData * getDecData() {return decData_;}

	int getNumSubproblems() {return decData_->getNumSubproblems();}

	int getNumCouplingRows() {return decData_->getNumCouplingRows();}

	int getNumCouplingCols() {return decData_->getNumCouplingCols();}

	int getNumSubproblemCouplingRows(int s) {return decData_->getNumSubproblemCouplingRows(s);}

	const int * getSubproblemCouplingRowIndices(int s) {return decData_->getSubproblemCouplingRowIndices(s);}

	int getNumSubproblemCouplingCols(int s) {return decData_->getNumSubproblemCouplingCols(s);}

	int getFullModelNumRows() {return det_->getNumRows() + decData_->getNumCouplingRows();}

	int getFullModelNumCols() {return det_->getNumCols();}

	int getNumIntegers() {return det_->getNumIntegers();}

	double evalLhsCouplingRow(int row, double ** solutions) {return decData_->evalLhsRow(row, solutions);}

	double evalLhsCouplingRowSubprob(int row, int subprob, double * subprobSolution) {return decData_->evalLhsRowSubprob(row, subprob, subprobSolution);}

	char getSenseCouplingRow(int row) {return decData_->getSenseRow(row);}

	double getRhsCouplingRow(int row) {return decData_->getRhsRow(row);}

	bool nonanticipativity() {return false;}

	bool isStochastic() {return false;}

	DSP_RTN_CODE decompose(
		int size,                     /**< [in] size of subproblem subset */
		int * subprobs,               /**< [in] subset of subproblems */
		int naux,                     /**< [in] number of auxiliary columns */
		double * clbd_aux,            /**< [in] lower bounds for auxiliary columns */
		double * cubd_aux,            /**< [in] upper bounds for auxiliary columns */
		double * obj_aux,             /**< [in] objective coefficients for auxiliary columns */
		CoinPackedMatrix *& mat,      /**< [out] constraint matrix */
		double *& clbd,               /**< [out] column lower bounds */
		double *& cubd,               /**< [out] column upper bounds */
		char   *& ctype,              /**< [out] column types */
		double *& obj,                /**< [out] objective coefficients */
		double *& rlbd,               /**< [out] row lower bounds */
		double *& rubd                /**< [out] row upper bounds */);

	DSP_RTN_CODE decomposeCoupling(
		int size,                    /**< [in] size of subproblem subset */
		int * subprobs,              /**< [in] subset of subproblems */
		CoinPackedMatrix *& cpl_mat, /**< [out] coupling constraint matrix */
		int *& cpl_cols,             /**< [out] columns of cpl_mat involved in coupling rows */
		int & cpl_ncols              /**< [out] size of cpl_cols */)
	{
		return decData_->decomposeCoupling(size, subprobs, cpl_mat, cpl_cols, cpl_ncols);
	}

	DSP_RTN_CODE getFullModel(
		CoinPackedMatrix *& mat, /**< [out] constraint matrix */
		double *& clbd,          /**< [out] column lower bounds */
		double *& cubd,          /**< [out] column upper bounds */
		char   *& ctype,         /**< [out] column types */
		double *& obj,           /**< [out] objective coefficients */
		double *& rlbd,          /**< [out] row lower bounds */
		double *& rubd           /**< [out] row upper bounds */);

	void __printData() {return det_->__printData();}

	/** only for setting decData_ lazily */
	void setDecData(DecData * decData)
	{
		assert(decData_ == NULL);
		decData_ = new DecData(*decData);
	}

	/** get objective coefficients */
	void getObjCoef(double * obj);

protected:

	DetModel * det_;       /**< underlying deterministic model */
	DecData * decData_;    /**< underlying decomposition information */

};

#endif /* DECDETMODEL_H_ */
