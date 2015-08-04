/*
 * DecModel.h
 *
 *  Created on: Jun 17, 2015
 *      Author: ctjandra
 */

#ifndef DECMODEL_H_
#define DECMODEL_H_

#include <map>
#include "CoinTime.hpp"
#include "SmiScnModel.hpp"

#include "Utility/StoMacros.h"
#include "Utility/StoRtnCodes.h"

/**
 * Interface for a decomposable model (not necessarily stochastic).
 * More specifically, this is a model with constraints of the form
 *    Ax + By = b, Hx = d,
 * where relaxing the latter decouples the problem into disjoint subproblems
 *    A^1 x^1 + B^1 y^1 = b^1, ..., A^k x^k + B^k y^k = b^k.
 * Here, Hx = d are the coupling constraints and x are the coupling variables.
 * Variables y and constraints that are not coupling are referred to as main.
 * (Currently coupling constraints must be equalities, but main constraints may be inequalities.)
 */
class DecModel
{
public:

	virtual ~DecModel() {}

public:

	/**
	 * Returns the number of subproblems this model decouples to.
	 */
	virtual int getNumSubproblems() = 0;

	/**
	 * Returns the number of coupling constraints (number of rows of H).
	 */
	virtual int getNumCouplingRows() = 0;

	/**
	 * Returns the number of coupling variables (number of columns of B or H, or dimension of x).
	 */
	virtual int getNumCouplingCols() = 0;

	/**
	 * Returns the number of coupling constraints for a subproblem s (number of rows of H involving x^s).
	 * When applying Lagrangian relaxation to Hx = d, this is the dimension of the Lagrangian multipliers
	 * restricted to the subproblem (i.e. multipliers from constraints irrelevant to the subproblem are ignored).
	 */
	virtual int getNumSubproblemCouplingRows(int s) = 0;

	/**
	 * Returns the indices of coupling constraints relevant to a subproblem s (indices of rows of H involving x^s).
	 * When applying Lagrangian relaxation to Hx = d, these are the indices of the Lagrangian multipliers
	 * restricted to the subproblem (i.e. multipliers from constraints irrelevant to the subproblem are ignored).
	 * (Note that these indices do not exist in the original problems; they are from 1 to getNumCouplingRows().)
	 */
	virtual const int * getSubproblemCouplingRowIndices(int s) = 0;

	/**
	 * Returns the number of coupling variables for a subproblem s (number of columns of A^s).
	 */
	virtual int getNumSubproblemCouplingCols(int s) = 0;

	/**
	 * Returns number of columns in full model. 
	 * (In stochastic, this is the extensive form. In deterministic, this is the complete model itself.)
	 */
	virtual int getFullModelNumCols() = 0;

	/**
	 * Returns the number of integer variables in the model.
	 */
	virtual int getNumIntegers() = 0;

	/**
	 * Stores in obj the objective coefficients of the model.
	 */
	virtual void getObjCoef(double * obj) = 0;

	/**
	 * Evaluate the left-hand side of a coupling row w.r.t. solutions in the subproblem space.
	 */
	virtual double evalLhsCouplingRow(int row, double ** solutions) = 0;

	/**
	 * Evaluate the left-hand side of a coupling row w.r.t. solution of one of the subproblems.
	 * Solution is in the space of the coupling variables of the subproblem, so this
	 * function is responsible for mapping the columns from the space of the subproblem
	 * to the space of the master problem.
	 */
	virtual double evalLhsCouplingRowSubprob(int row, int subprob, double * subprobSolution) = 0;

	/**
	 * Return the sense of a coupling row.
	 */
	virtual char getSenseCouplingRow(int row) = 0;

	/**
	 * Return the right-hand side of a coupling row.
	 */
	virtual double getRhsCouplingRow(int row) = 0;

	/**
	 * Indicates whether the coupling constraints are nonanticipativity constraints
	 * (i.e. x_1 = ... = x_s), which may receive special treatment due to their structure.
	 * If nonanticipativity is true, implementation must ensure all these are equal:
	 *   the dimension of x_s for any s,
	 *   getNumCouplingRows() / getNumSubproblems(),
	 *   getNumCouplingCols() / getNumSubproblems(),
	 *   getNumSubproblemCouplingRows(s) for any s, and
	 *   getNumSubproblemCouplingCols(s) for any s.
	 * See class DecTssModel for more details on the special treatment for nonanticipativity.
	 */
	virtual bool nonanticipativity() = 0;

	/**
	 * If true, this is a stochastic model that can be downcasted to StoModel.
	 * Used to handle specific stochastic cases.
	 */
	virtual bool isStochastic() = 0;

	/**
	 * This routine decomposes the model based on inputs. If size = 0, then
	 * this results in a standard Benders decomposition structure. If size = 1, then
	 * this results in a dual decomposition structure.
	 */
	virtual STO_RTN_CODE decompose(
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
		double *& rubd               /**< [out] row upper bounds */) = 0;

	/**
	 * Produces the coupling constraints for the given subproblems. The coupling constraint matrix
	 * includes only what will be sent to a subproblem. More precisely, it contains the coupling
	 * constraints that involve at least one variable from a subproblem in subprobs. Its column
	 * dimension must be the sum of getNumSubproblemCouplingCols(s) for all s in subprobs. The
	 * array cpl_cols is the subset of columns of cpl_mat that are involved in coupling rows,
	 * which will be sent from a subproblem to the master).
	 */
	virtual STO_RTN_CODE decomposeCoupling(
		int size,                    /**< [in] size of subproblem subset */
		int * subprobs,              /**< [in] subset of subproblems */
		CoinPackedMatrix *& cpl_mat, /**< [out] coupling constraint matrix */
		int *& cpl_cols,             /**< [out] columns of cpl_mat involved in coupling rows */
		int & cpl_ncols              /**< [out] size of cpl_cols */) = 0;

	/**
	 * Returns the full model in matrix form, including coupling constraints.
	 */
	virtual STO_RTN_CODE getFullModel(
		CoinPackedMatrix *& mat, /**< [out] constraint matrix */
		double *& clbd,          /**< [out] column lower bounds */
		double *& cubd,          /**< [out] column upper bounds */
		char   *& ctype,         /**< [out] column types */
		double *& obj,           /**< [out] objective coefficients */
		double *& rlbd,          /**< [out] row lower bounds */
		double *& rubd           /**< [out] row upper bounds */) = 0;

	/**
	 * Print model nicely.
	 */
	virtual void __printData() = 0;
};

#endif /* DECMODEL_H_ */
