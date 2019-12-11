/*
 * DecTssModel.h
 *
 *  Created on: Jun 17, 2015
 *      Author: ctjandra
 */

#ifndef DECTSSMODEL_H_
#define DECTSSMODEL_H_

#include <numeric>
#include "Model/TssModel.h"
#include "Model/DecModel.h"

/**
 * Decomposable version of the two-stage stochastic problem through nonanticipativity constraints.
 *
 * Nonanticipativity constraints are represented as x = x_i for every scenario i.
 * The variable x is eliminated from the Lagrangian relaxation by adding the constraint that the
 * sum of the Lagrange multipliers over scenarios is zero (see e.g. [1] for more precise details).
 * The master problem should either:
 * 1. take this into account by adding the normalization constraint, or
 * 2. use the methods to convert to the alternative representation x_i = x_{i+1} for all i.
 *
 * [1] M Lubin, K Martin, CG Petra, B Sandıkçı, "On parallelizing dual decomposition in stochastic 
 *     integer programming", 2013
 */
class DecTssModel: public TssModel, public DecModel {

public:

	/** default constructor */
	DecTssModel();

	/** copy constructor */
	DecTssModel(const DecTssModel & rhs);

	/** copy constructor */
	DecTssModel(const TssModel & rhs);

	/** default destructor */
	virtual ~DecTssModel();

public:

	/**
	 * Returns the set of subproblem indices for a given coupling column j.
	 */
	std::vector<int> getCoupledSubproblemIndices(int j) {
		std::vector<int> cols(getNumSubproblems());
		std::iota(cols.begin(), cols.begin() + getNumSubproblems(), 0);
		return cols;
	}
	/**
	 * Returns the number of scenarios.
	 */
	int getNumSubproblems() {return getNumScenarios();}

	/**
	 * Returns the number of first-stage variables.
	 */
	int getNumCouplingCols() {return getNumCols(0);}

	/**
	 * Returns the number of first-stage variables.
	 */
	int getNumSubproblemCouplingCols(int s) {return getNumCols(0);}

	/**
	 * Returns the first-stage variable indices.
	 */
	const int * getSubproblemCouplingColIndices(int s);

	/**
	 * Returns the number of columns in the extensive form.
	 */
	int getFullModelNumRows() {return getNumRows(0) + getNumScenarios() * getNumRows(1);}

	/**
	 * Returns the number of columns in the extensive form.
	 */
	int getFullModelNumCols() {return getNumCols(0) + getNumScenarios() * getNumCols(1);}

	/**
	 * Returns the number of integer variables in the extensive form.
	 */
	int getNumIntegers() {return TssModel::getNumIntegers(0) + getNumScenarios() * TssModel::getNumIntegers(1);}

	int getNumCouplingIntegers() {return TssModel::getNumIntegers(0) * getNumScenarios();}

	double evalLhsCouplingRow(int row, double ** solutions);

	/**@name Functions specific for the Dual Decomposition */
	//@{

	/**
	 * This considers a dual decomposition and returns the number of non-anticipativity constraints.
	 */
	int getNumCouplingRows() {return getNumScenarios() * getNumCols(0);}

	/**
	 * This considers a dual decomposition and returns the number of first-stage variables copied for each scenario s.
	 */
	int getNumSubproblemCouplingRows(int s) {return getNumCols(0);}

	/**
	 * This returns the row-th element of vector lambda_j^T x_j.
	 * @param row is in [0, getNumCouplingRows()].
	 * @param subprob is a scenario index (j) ranged in [0, getNumSubproblems()].
	 * @param subprobSolution is the first-stage solution (x) of the size getNumCouplingCols().
	 */
	double evalLhsCouplingRowSubprob(int row, int subprob, double * subprobSolution);

	/**
	 * Returns the coupling row lower bound.
	 */
	double getCouplingRowLower(int row) {return 0;}

	/**
	 * Returns the coupling row upper bound.
	 */
	double getCouplingRowUpper(int row) {return 0;}

	/**
	 * Returns the coupling row sense.
	 */
	char getSenseCouplingRow(int row) {return 'E';}

	/**
	 * Returns the coupling row rhs.
	 */
	double getRhsCouplingRow(int row) {return 0;}

	//@}

	bool nonanticipativity() {return true;}

	bool isStochastic() {return true;}

	DSP_RTN_CODE decompose(
		int size,                    /**< [in] size of subproblem subset */
		int * subprobs,              /**< [in] subset of subproblems */
		int naux,                    /**< [in] number of auxiliary columns */
		double * clbd_aux,           /**< [in] lower bounds for auxiliary columns */
		double * cubd_aux,           /**< [in] upper bounds for auxiliary columns */
		double * obj_aux,            /**< [in] objective coefficients for auxiliary columns */
		CoinPackedMatrix *& mat,     /**< [out] constraint matrix */
		double *& clbd,              /**< [out] column lower bounds */
		double *& cubd,              /**< [out] column upper bounds */
		char   *& ctype,             /**< [out] column types */
		double *& obj,               /**< [out] objective coefficients */
		double *& rlbd,              /**< [out] row lower bounds */
		double *& rubd               /**< [out] row upper bounds */)
	{
		return TssModel::decompose(size, subprobs, naux, clbd_aux, cubd_aux, obj_aux, mat,
			clbd, cubd, ctype, obj, rlbd, rubd);
	}

	DSP_RTN_CODE decomposeCoupling(
		int size,                    /**< [in] size of subproblem subset */
		int * subprobs,              /**< [in] subset of subproblems */
		CoinPackedMatrix *& cpl_mat, /**< [out] coupling constraint matrix */
		int *& cpl_cols,             /**< [out] columns of cpl_mat involved in coupling rows */
		int & cpl_ncols              /**< [out] size of cpl_cols */);

	/**
	 * This creates the subproblems with coupling columns; that is,
	 *   lb^k <= A^k x^k + B^k y^k <= ub^k
	 */
	DSP_RTN_CODE copySubprob(
		int subprob,             /**< [in] subproblem index */
		CoinPackedMatrix *& mat, /**< [out] constraint matrix [A^k B^k] */
		double *& clbd,          /**< [out] column lower bounds of y */
		double *& cubd,          /**< [out] column upper bounds of y */
		char   *& ctype,         /**< [out] column types of y */
		double *& obj,           /**< [out] objective coefficients for y */
		double *& rlbd,          /**< [out] row lower bounds */
		double *& rubd           /**< [out] row upper bounds */);

	DSP_RTN_CODE copyRecoProb(
		int scen,                     /**< [in] scenario index */
		CoinPackedMatrix *& mat_tech, /**< [out] technology matrix (A matrix) */
		CoinPackedMatrix *& mat_reco, /**< [out] recourse matrix (B matrix) */
		double *& clbd_reco,          /**< [out] column lower bounds of y */
		double *& cubd_reco,          /**< [out] column upper bounds of y */
		char   *& ctype_reco,         /**< [out] column types of y */
		double *& obj_reco,           /**< [out] objective coefficients for y */
		double *& rlbd_reco,          /**< [out] row lower bounds */
		double *& rubd_reco           /**< [out] row upper bounds */);

	DSP_RTN_CODE getFullModel(
		CoinPackedMatrix *& mat, /**< [out] constraint matrix */
		double *& clbd,          /**< [out] column lower bounds */
		double *& cubd,          /**< [out] column upper bounds */
		char   *& ctype,         /**< [out] column types */
		double *& obj,           /**< [out] objective coefficients */
		double *& rlbd,          /**< [out] row lower bounds */
		double *& rubd           /**< [out] row upper bounds */);

	/** Methods on alternative representation used in subgradient algorithm */

	/** Evaluates row with nonanticipativity represented as x_i = x_{i+1} */
	double evalLhsCouplingRowAlternative(int row, double ** solutions);

	/** Converts Lagrangian multipliers back from the alternative representation x_i = x_{i+1} (to x = x_i) */
	void convertLagrangianFromAlternative(double * multipliers, double *& newMultipliers);

	void __printData() {return TssModel::__printData();}

private:

	int* master_col_indices_; /**< master column indices */
};

#endif /* DECTSSMODEL_H_ */
