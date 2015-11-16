/*
 * DecTssModel.h
 *
 *  Created on: Jun 17, 2015
 *      Author: ctjandra
 */

#ifndef DECTSSMODEL_H_
#define DECTSSMODEL_H_

#include "TssModel.h"
#include "DecModel.h"

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

	int getNumSubproblems() {return getNumScenarios();}

	int getNumCouplingRows() {
		return getNumScenarios() * getNumCols(0); /* rows of nonanticipativity constraints */
	}

	int getNumCouplingCols() {
		return getNumScenarios() * getNumCols(0); /* copies of first-stage variables for all scenarios */
	}

	int getNumSubproblemCouplingRows(int s) {
		return getNumCols(0); /* copies of first-stage variables for scenario s */
	}

	const int * getSubproblemCouplingRowIndices(int s);

	int getNumSubproblemCouplingCols(int s) {return getNumCols(0);}

	int getFullModelNumCols() {return getNumCols(0) + getNumScenarios() * getNumCols(1);}

	int getNumIntegers() {return getNumCoreIntegers();}

	double evalLhsCouplingRow(int row, double ** solutions);

	double evalLhsCouplingRowSubprob(int row, int subprob, double * subprobSolution);

	char getSenseCouplingRow(int row) {return 'E';}

	double getRhsCouplingRow(int row) {return 0;}

	bool nonanticipativity() {return true;}

	bool isStochastic() {return true;}

	STO_RTN_CODE decompose(
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

	STO_RTN_CODE decomposeCoupling(
		int size,                    /**< [in] size of subproblem subset */
		int * subprobs,              /**< [in] subset of subproblems */
		CoinPackedMatrix *& cpl_mat, /**< [out] coupling constraint matrix */
		int *& cpl_cols,             /**< [out] columns of cpl_mat involved in coupling rows */
		int & cpl_ncols              /**< [out] size of cpl_cols */);

	STO_RTN_CODE getFullModel(
		CoinPackedMatrix *& mat, /**< [out] constraint matrix */
		double *& clbd,          /**< [out] column lower bounds */
		double *& cubd,          /**< [out] column upper bounds */
		char   *& ctype,         /**< [out] column types */
		double *& obj,           /**< [out] objective coefficients */
		double *& rlbd,          /**< [out] row lower bounds */
		double *& rubd           /**< [out] row upper bounds */);

	/** get objective coefficients */
	void getObjCoef(double * obj);


	/** Methods on alternative representation used in subgradient algorithm */

	/** Evaluates row with nonanticipativity represented as x_i = x_{i+1} */
	double evalLhsCouplingRowAlternative(int row, double ** solutions);

	/** Converts Lagrangian multipliers back from the alternative representation x_i = x_{i+1} (to x = x_i) */
	void convertLagrangianFromAlternative(double * multipliers, double *& newMultipliers);

	void __printData() {return TssModel::__printData();}
};

#endif /* DECTSSMODEL_H_ */
