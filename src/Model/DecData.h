/*
 * DecData.h
 *
 *  Created on: Jun 17, 2015
 *      Author: ctjandra
 */

#ifndef DECDATA_H_
#define DECDATA_H_

#include "TssModel.h"
#include "DecModel.h"

/** 
 * Data structure to explicitly store coupling constraints and column partition, meant to be attached
 * to a decomposition model.
 */
class DecData {

public:

	/** default constructor */
	DecData(
		int      nsubprobs,      /**< number of subproblems */
		int      ncols,          /**< number of columns for the entire problem */
		int      ncoupling,      /**< number of coupling constraints */
		int *    varPartition,   /**< partition of columns into subproblems */
		int *    couplingStarts, /**< indices in cols at which each coupling constraint starts */
		int *    couplingCols,   /**< variables of each coupling constraint left-hand side */
		double * couplingCoeffs, /**< coefficients of each coupling constraint left-hand side */
		char *   couplingSenses, /**< senses of each coupling constraint */
		double * couplingRhs     /**< right-hand sides of each coupling constraint */);

	/** copy constructor */
	DecData(const DecData & rhs);

	/** default destructor */
	virtual ~DecData();

public:

	int getNumSubproblems() {return nsubprobs_;}

	int getNumCouplingRows() {return nrows_coupling_;}

	int getNumCouplingCols() {return ncols_coupling_;}

	int getNumCols() {return ncols_;}

	int getNumSubproblemCouplingRows(int s) {return nrows_subprob_coupling_[s];}

	const int * getSubproblemCouplingRowIndices(int s) {return rows_subprob_[s];}

	int getNumSubproblemCouplingCols(int s) {return ncols_subprob_coupling_[s];}

	const int * getColPartition() {return cols_subprob_;}

	const CoinPackedMatrix * getConstraintMatrix() {return rows_coupling_;}

	/** Evaluate the left-hand side of a row w.r.t. solutions in the subproblem space */
	double evalLhsRow(int row, double ** solutions);

	/** Evaluate the left-hand side of a row w.r.t. solution of one of the subproblems */
	double evalLhsRowSubprob(int row, int subprob, double * subprobSolution);

	/** Return the sense of a row */
	char getSenseRow(int row) {return senses_[row];}

	/** Return the right-hand side of a row */
	double getRhsRow(int row) {return rhs_coupling_[row];}

	/**
	 * Helper function to decompose a problem given in the form of a matrix.
	 * Given a matrix from the entire problem, selects columns corresponding to the given subproblems
	 * and discards rows that correspond to other subproblems (i.e. not involving any selected columns).
	 */
	STO_RTN_CODE decompose(
		int size,                /**< [in] size of subproblem subset */
		int * subprobs,          /**< [in] subset of subproblems */
		CoinPackedMatrix *& mat, /**< [in/out] constraint matrix */
		double *& clbd,          /**< [in/out] column lower bounds */
		double *& cubd,          /**< [in/out] column upper bounds */
		char   *& ctype,         /**< [in/out] column types */
		double *& obj,           /**< [in/out] objective coefficients */
		double *& rlbd,          /**< [in/out] row lower bounds */
		double *& rubd           /**< [in/out] row upper bounds */);

	/** Helper function to decompose the coupling constraint matrix into subproblems */
	STO_RTN_CODE decomposeCoupling(
		int size,                    /**< [in] size of subproblem subset */
		int * subprobs,              /**< [in] subset of subproblems */
		CoinPackedMatrix *& cpl_mat, /**< [out] coupling constraint matrix */
		int *& cpl_cols,             /**< [out] columns of cpl_mat involved in coupling rows */
		int & cpl_ncols              /**< [out] size of cpl_cols */);

private:

	/** Check if the subproblem of a given column is in a list of subproblems */
	bool columnInSubprobs(
		int col,                 /**< column to check if in subproblem subset */
		int size,                /**< size of subproblem subset */
		int * subprobs           /**< subset of subproblems */);

	/** Select the columns of a matrix corresponding to the given subproblems and eliminate zero rows */
	void reduceMatrixToSubproblems(
		int size,                    /**< [in] size of subproblem subset */
		int * subprobs,              /**< [in] subset of subproblems */
		CoinPackedMatrix *& mat,     /**< [in/out] constraint matrix */
		vector<int> & rowsDiscarded, /**< [out] row indices from the matrix that were discarded */
		vector<int> & colsDiscarded  /**< [out] column indices from the matrix that were discarded */);

	/**
	 * Return false if there are two nonzero coefficients of different classes in the same row
	 * of a given matrix.
	 */
	bool checkPartitionConsistency(const CoinPackedMatrix * mat);

protected:

	int nsubprobs_;           /**< number of subproblems */
	int nrows_coupling_;      /**< number of coupling constraints */
	int ncols_coupling_;      /**< number of coupling variables */
	int * nrows_subprob_coupling_; /**< number of coupling constraints for each subproblem */
	int * ncols_subprob_coupling_; /**< number of coupling variables for each subproblem */
	int * ncols_subprob_;     /**< number of variables for each subproblem */

	int ncols_;               /**< number of columns of the entire problem */
	int * cols_subprob_;      /**< partition of columns of the entire problem into subproblems */
	int ** rows_subprob_;     /**< lists of row indices of rows_coupling_ for each subproblem */
	bool * cols_incoupling_;  /**< indicates whether a column is part of coupling constraints */

	int ** orig_coupling_col_; /**< orig_coupling_col_[s][i] maps the i-th coupling variable of 
	                             *  subproblem s to the variable index in the original space */

	CoinPackedMatrix * rows_coupling_; /**< coupling constraints (on original space of variables) */
	char * senses_;            /** row senses for coupling constraints ([L]ess than or equal, [E]qual, or [G]reater than or equal) */
	double * rhs_coupling_;    /** right-hand side of coupling constraints (on the original space of variables) */


};

#endif /* DECDATA_H_ */
