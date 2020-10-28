/*
 * DecModel.h
 *
 *  Created on: Jun 17, 2015
 *      Author: ctjandra
 */

#ifndef DECMODEL_H_
#define DECMODEL_H_

#include <map>
#include <vector>
/** Coin */
#include "CoinTime.hpp"
#include "CoinPackedMatrix.hpp"
/** Dsp */
#include "Utility/DspMacros.h"
#include "Utility/DspRtnCodes.h"


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

/**
 * Extension:
 * The design follows the matrix of the form
 *   lb^0 <= H^1 x^1 + ... + H^k x^k                     <= ub^0
 *   lb^1 <= A^1 x^1                 + B^1 y^1           <= ub^1
 *   ...
 *   lb^k <=                 A^k x^k           + B^k y^k <= ub^k
 *
 * The nonzero elements of the constraint matrix look like this:
 *   X X X       <- coupling rows
 *   X     X
 *     X     X
 *       X     X
 *   ^ ^ ^
 *   coupling columns
 *
 * NOTE: A Benders-type method will fix the coupling column values,
 *       whereas a Lagrangian method will relax the coupling rows.
 *
 * Either fixing the coupling column values or relaxing the coupling rows
 * decomposes the problem into three subproblems.
 */
class DecModel
{
public:

	virtual ~DecModel() {}

public:

	/**
	 * Returns the set of subproblem indices for a given coupling column j.
	 */
	virtual std::vector<int> getCoupledSubproblemIndices(int j) = 0;

	/**
	 * Returns the number of subproblems this model decouples to.
	 */
	virtual int getNumSubproblems() = 0;

	/**
	 * Returns the number of coupling constraints (number of rows of H).
	 */
	virtual int getNumCouplingRows() = 0;

	/**
	 * Returns the number of coupling variables (dimension of [x^1, ..., x^k]).
	 */
	virtual int getNumCouplingCols() = 0;

	/**
	 * Returns the number of coupling constraints for a subproblem s (number of rows of H involving x^s).
	 * When applying Lagrangian relaxation to Hx = d, this is the dimension of the Lagrangian multipliers
	 * restricted to the subproblem (i.e. multipliers from constraints irrelevant to the subproblem are ignored).
	 */
	virtual int getNumSubproblemCouplingRows(int s) = 0;

	/**
	 * Returns the number of coupling columns for subproblem s (number of columns of H^s)
	 */
	virtual int getNumSubproblemCouplingCols(int s) = 0;

	/**
	 * Returns the indices of coupling columns for subproblem s (column indices of H^s)
	 * with respect to the full problem
	 */
	virtual const int * getSubproblemCouplingColIndices(int s) = 0;

	/**
	 * Returns number of rows in full model.
	 * (In stochastic, this is the extensive form. In deterministic, this is the complete model itself.)
	 */
	virtual int getFullModelNumRows() = 0;

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
	 * Returns the number of integer coupling variables
	 */
	virtual int getNumCouplingIntegers() = 0;

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
	 * Retruns the coupling column objective coefficients
	 */
	virtual const double * getCouplingColsObjs() = 0;

	/**
	 * Returns the coupling row lower bound
	 */
	virtual double getCouplingRowLower(int row) = 0;

	/**
	 * Returns the coupling row upper bound
	 */
	virtual double getCouplingRowUpper(int row) = 0;

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
	 * If true, this models contains quadratic terms in the objectives or constraints
	 **/
	virtual bool isQCQP() = 0;
	
	/**
	 * Set whether this is dro model or not.
	 */
	virtual void setDro(bool yes) = 0;

	/**
	 * If true, this is a distributionally robust variant.
	 */
	virtual bool isDro() = 0;

	/**
	 * Returns the number of reference scenarios (for DRO).
	 */
	virtual int getNumReferences() = 0;

	/**
	 * Returns the Wasserstein distance limit for DRO.
	 */
	virtual double getWassersteinSize() = 0;

	/**
	 * Returns the Wasserstein distance between two random realizations i and j.
	 */
	virtual double getWassersteinDist(int i, int j) = 0;

	/**
	 * Returns the probability of reference i.
	 */
	virtual double getReferenceProbability(int i) = 0;

	/**
	 * This routine decomposes the model based on inputs. If size = 0, then
	 * this results in a standard Benders decomposition master problem. If size = 1, then
	 * this results in a dual decomposition subproblem.
	 */
	virtual DSP_RTN_CODE decompose(
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

	virtual DSP_RTN_CODE decompose(
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
		CoinPackedMatrix *& qobj,     /**< [out] quadratic objective */
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
	virtual DSP_RTN_CODE decomposeCoupling(
		int size,                    /**< [in] size of subproblem subset */
		int * subprobs,              /**< [in] subset of subproblems */
		CoinPackedMatrix *& cpl_mat, /**< [out] coupling constraint matrix */
		int *& cpl_cols,             /**< [out] columns of cpl_mat involved in coupling rows */
		int & cpl_ncols              /**< [out] size of cpl_cols */) = 0;

	/**
	 * This creates the subproblems with coupling columns; that is,
	 *   lb^k <= A^k x^k + B^k y^k <= ub^k
	 */
	virtual DSP_RTN_CODE copySubprob(
		int subprob,             /**< [in] subproblem index */
		CoinPackedMatrix *& mat, /**< [out] constraint matrix [A^k B^k] */
		double *& clbd,          /**< [out] column lower bounds of y */
		double *& cubd,          /**< [out] column upper bounds of y */
		char   *& ctype,         /**< [out] column types of y */
		double *& obj,           /**< [out] objective coefficients for y */
		double *& rlbd,          /**< [out] row lower bounds */
		double *& rubd           /**< [out] row upper bounds */) = 0;

	/**
	 * This creates recourse problem structure for a given scenario index; that is,
	 *   lb^k <= B^k y^k <= ub^k with A^k matrix separately.
	 */
	virtual DSP_RTN_CODE copyRecoProb(
		int scen,                      /**< [in] scenario index */
		CoinPackedMatrix *& mat_tech,  /**< [out] technology matrix (A matrix) */
		CoinPackedMatrix *& mat_reco,  /**< [out] recourse matrix (B matrix) */
		double *& clbd_reco,           /**< [out] column lower bounds of y */
		double *& cubd_reco,           /**< [out] column upper bounds of y */
		char   *& ctype_reco,          /**< [out] column types of y */
		double *& obj_reco,            /**< [out] objective coefficients for y */
		double *& rlbd_reco,           /**< [out] row lower bounds */
		double *& rubd_reco,           /**< [out] row upper bounds */
		bool adjust_probability = true /**< [in] adjust probability (only for stochastic)*/) = 0;

		virtual DSP_RTN_CODE copyRecoProb(
		int scen,                     /**< [in] scenario index */
		CoinPackedMatrix *& mat_tech, /**< [out] technology matrix (A matrix) */
		CoinPackedMatrix *& mat_reco, /**< [out] recourse matrix (B matrix) */
		double *& clbd_reco,          /**< [out] column lower bounds of y */
		double *& cubd_reco,          /**< [out] column upper bounds of y */
		char   *& ctype_reco,         /**< [out] column types of y */
		double *& obj_reco,           /**< [out] objective coefficients for y */
		CoinPackedMatrix *& qobj_reco_coupling,/**< [out] coupling quadratric coefficients (y^2}*/
		CoinPackedMatrix *& qobj_reco_ncoupling, /**< [out] non-coupling quadratic coefficients (xy) */
		double *& rlbd_reco,          /**< [out] row lower bounds */
		double *& rubd_reco           /**< [out] row upper bounds */) = 0;

	/**
	 * Returns the full model in matrix form, including coupling constraints.
	 */
	virtual DSP_RTN_CODE getFullModel(
		CoinPackedMatrix *& mat, /**< [out] constraint matrix */
		double *& clbd,          /**< [out] column lower bounds */
		double *& cubd,          /**< [out] column upper bounds */
		char   *& ctype,         /**< [out] column types */
		double *& obj,           /**< [out] objective coefficients */
		double *& rlbd,          /**< [out] row lower bounds */
		double *& rubd           /**< [out] row upper bounds */) = 0;
	virtual DSP_RTN_CODE getFullModel(
		CoinPackedMatrix *& mat, /**< [out] constraint matrix */
		double *& clbd,          /**< [out] column lower bounds */
		double *& cubd,          /**< [out] column upper bounds */
		char   *& ctype,         /**< [out] column types */
		double *& obj,           /**< [out] objective coefficients */
		CoinPackedMatrix *& qobj,/**< [out] quadratic objective coefficients */
		double *& rlbd,          /**< [out] row lower bounds */
		double *& rubd           /**< [out] row upper bounds */) = 0;
	/**
	 * Print model nicely.
	 */
	virtual void __printData() = 0;
};

#endif /* DECMODEL_H_ */
