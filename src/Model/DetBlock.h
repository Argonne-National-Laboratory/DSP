/*
 * DetBlock.h
 *
 *  Created on: Aug 31, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_MODEL_DETBLOCK_H_
#define SRC_MODEL_DETBLOCK_H_

#include "DetModel.h"

class DetBlock: public DetModel {
public:

	/** default constructor */
	DetBlock() :
		DetModel(),
		num_coupling_cols_(0),
		num_coupling_rows_(0),
		coupling_cols_(NULL),
		coupling_rows_(NULL) {
		/** nothing to do */
	}

	DetBlock(
			const CoinPackedMatrix * mat, /**< constraint matrix */
			const double * clbd,          /**< column lower bounds */
			const double * cubd,          /**< column upper bounds */
			const char   * ctype,         /**< column types */
			const double * obj,           /**< objective coefficients */
			const double * rlbd,          /**< row lower bounds */
			const double * rubd           /**< row upper bounds */) :
				DetModel(mat,clbd,cubd,ctype,obj,rlbd,rubd),
				num_coupling_cols_(0),
				num_coupling_rows_(0),
				coupling_cols_(NULL),
				coupling_rows_(NULL) {
		/** nothing to do */
	}

	DetBlock(
			const CoinPackedMatrix * mat, /**< constraint matrix */
			const double * clbd,          /**< column lower bounds */
			const double * cubd,          /**< column upper bounds */
			const char   * ctype,         /**< column types */
			const double * obj,           /**< objective coefficients */
			const CoinPackedMatrix * qobj,/**< quadratic objective coefficients */
			const double * rlbd,          /**< row lower bounds */
			const double * rubd           /**< row upper bounds */) :
				DetModel(mat,clbd,cubd,ctype,obj,qobj,rlbd,rubd),
				num_coupling_cols_(0),
				num_coupling_rows_(0),
				coupling_cols_(NULL),
				coupling_rows_(NULL) {
		/** nothing to do */
	}

	DetBlock(
			const CoinBigIndex * start, /**< start index for each row */
			const int    * index,       /**< column indices */
			const double * value,       /**< constraint elements */
			const int      numels,      /**< number of elements in index and value */
			const int      ncols,       /**< number of columns */
			const int      nrows,       /**< number of rows */
			const double * clbd,        /**< column lower bounds */
			const double * cubd,        /**< column upper bounds */
			const char   * ctype,       /**< column types */
			const double * obj,         /**< objective coefficients */
			const double * rlbd,        /**< row lower bounds */
			const double * rubd         /**< row upper bounds */) :
				DetModel(start,index,value,numels,ncols,nrows,clbd,cubd,ctype,obj,rlbd,rubd),
				num_coupling_cols_(0),
				num_coupling_rows_(0),
				coupling_cols_(NULL),
				coupling_rows_(NULL) {
		/** nothing to do */
	}

	DetBlock(
			const CoinBigIndex * start, /**< start index for each row */
			const int    * index,       /**< column indices */
			const double * value,       /**< constraint elements */
			const int      numels,      /**< number of elements in index and value */
			const int      ncols,       /**< number of columns */
			const int      nrows,       /**< number of rows */
			const double * clbd,        /**< column lower bounds */
			const double * cubd,        /**< column upper bounds */
			const char   * ctype,       /**< column types */
			const double * obj,         /**< objective coefficients */
			const CoinBigIndex * qobjstart, /**< quadratic objective start index for each row */
			const int	 * qobjindex,	/**< quadratic objective index */
			const double * qobjvalue,	/**< quadratic objective coefficient*/
			const int	   qobjnumels,	/**< number of elements in qobj index and value */
			const double * rlbd,        /**< row lower bounds */
			const double * rubd         /**< row upper bounds */) :
				DetModel(start,index,value,numels,ncols,nrows,clbd,cubd,ctype,obj,qobjstart,qobjindex,qobjvalue,qobjnumels,rlbd,rubd),
				num_coupling_cols_(0),
				num_coupling_rows_(0),
				coupling_cols_(NULL),
				coupling_rows_(NULL) {
		/** nothing to do */
	}

	/** copy constructor */
	DetBlock(const DetBlock& rhs);

	/** default destructor */
	virtual ~DetBlock() {
		FREE_ARRAY_PTR(coupling_cols_);
		FREE_ARRAY_PTR(coupling_rows_);
	}

	/** get number of coupling columns */
	int getNumCouplingCols() {return num_coupling_cols_;}

	/** get number of coupling rows */
	int getNumCouplingRows() {return num_coupling_rows_;}

	/** get coupling columns */
	int* getCouplingCols() {return coupling_cols_;}

	/** get coupling rows */
	int* getCouplingRows() {return coupling_rows_;}

	/** set coupling columns */
	void setCouplingCols(int n, int* index) {
		FREE_ARRAY_PTR(coupling_cols_);
		coupling_cols_ = new int [n];
		CoinCopyN(index, n, coupling_cols_);
		num_coupling_cols_ = n;
	}

	/** set coupling columns */
	void setCouplingRows(int n, int* index) {
		FREE_ARRAY_PTR(coupling_rows_);
		coupling_rows_ = new int [n];
		CoinCopyN(index, n, coupling_rows_);
		num_coupling_rows_ = n;
	}

protected:

	int num_coupling_cols_; /**< number of coupling columns */
	int num_coupling_rows_; /**< number of coupling rows */
	int* coupling_cols_;    /**< coupling columns */
	int* coupling_rows_;    /**< coupling rows */
};

#endif /* SRC_MODEL_DETBLOCK_H_ */
