/*
 * DdMasterDsb.h
 *
 *  Created on: Feb 16, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DUALDECOMP_DDMASTERDSB_H_
#define SRC_SOLVER_DUALDECOMP_DDMASTERDSB_H_

#include "Solver/DualDecomp/DdMaster.h"

/** Implementation of Doubly Stabilized Bundle Method
 *
 * General problem description:
 *   min -l^k v + \sum_s \sum_j [ d_s(\lambda^j) + (B_s x_s^j)^T (prox^k - \lambda^j) ] y_s^j
 *       + 0.5 tau_k \sum_s | \sum_j B_s x_s^j y_s^j + b u |_2^2
 *   subject to
 *     u - \sum_j y_s^j = 0 \forall s
 *     u - v = 1
 *     b u + \sum_j \sum_s B_s x_s^j y_s^j \geq -prox^k / \tau_k (**)
 *     u. v, y_s^j >= 0,
 *
 * where constraint (**) is considered for the inequality coupling constraint
 *   \sum_s B_s x_s \leq b
 * and the inequality sign of (**) should depend on that of the coupling constraint.
 *
 * Lower triangle of the Hessian is:
 *   \tau_k * [ |S||b|^2      ,                           ,                      ,              ,               ; (u)
 *               0.0          , 0.0,                      ,                      ,              ,               ; (v)
 *               b^T B_1 x_1^1, 0.0, |B_1 x_1^1|^2        ,                      ,              ,               ; (y_1^1)
 *               b^T B_2 x_2^1, 0.0, 0.0                  , |B_1 x_2^1|^2        ,              ,               ; (y_2^1)
 *               b^T B_1 x_1^2, 0.0, x_1^1 B_1^T B_1 x_1^2, 0.0                  , |B_1 x_1^2|^2,               ; (y_1^2)
 *               b^T B_2 x_2^2, 0.0, 0.0                  , x_2^1 B_2^T B_2 x_2^2, 0.0          , |B_2 x_2^2|^2 ] (y_2^2)
 *
 *
 * Stochastic problem description (decoupling all Lagrangian multipliers):
 *   min -l^k v + \sum_s \sum_j [ d_s(\lambda_s^j) + (x_s^j)^T (prox_s^k - \lambda_s^j) ] y_s^j
 *       + 0.5 tau_k \sum_s | \sum_j x_s^j y_s^j + w |_2^2
 *   subject to
 *     u - \sum_j y_s^j = 0 \forall s
 *     u - v = 1
 *     u, v, y_s^j >= 0, w \in R^n
 *
 * Lower triangle of the Hessian is:
 *   \tau_k * [ 0.0,                            ,                 ,           ,           ; (u)
 *              0.0, 0.0,                       ,                 ,           ,           ; (v)
 *              0.0, 0.0, |S|  ,                ,                 ,           ,           ; (w)
 *              0.0, 0.0, x_1^1, |x_1^1|^2      ,                 ,           ,           ; (y_1^1)
 *              0.0, 0.0, x_2^1, 0.0            , |x_2^1|^2       ,           ,           ; (y_2^1)
 *              0.0, 0.0, x_1^2, (x_1^2)^T x_1^1, 0.0             , |x_1^2|^2 ,           ; (y_1^2)
 *              0.0, 0.0, x_2^2, 0.0            , (x_2^2)^T x_2^1 , 0.0       , |x_2^2|^2 ] (y_2^2)
 */
class DdMasterDsb: public DdMaster {
public:

	/** constructor */
	DdMasterDsb(
			DspParams *  par,    /**< parameter pointer */
			DecModel *   model,  /**< model pointer */
			DspMessage * message /**< message pointer */);

	/** destructor */
	virtual ~DdMasterDsb();

	/** initialize */
	virtual DSP_RTN_CODE init();

	/** solve */
	virtual DSP_RTN_CODE solve();

	/** update problem */
	virtual DSP_RTN_CODE updateProblem();

protected:

	/** create problem */
	virtual DSP_RTN_CODE createProblem();

private:

	/** manage bundle */
	DSP_RTN_CODE manageBundle(
			double * objvals /**< subproblem objective values */);

	/** apply level change */
	DSP_RTN_CODE applyLevelChange();

	/** update dynamic objective coefficients */
	DSP_RTN_CODE updateDynObjs();

private:

	double * prox_;      /**< proximal point */
	double * lambda_;    /**< dual variable */
	double phi_t_;       /**< predicted increase from bundle */
	double phi_l_;       /**< target increase from level */
	double alpha_t_;     /**< bundle test parameter */
	double alpha_l_;     /**< level update parameter */
	double eps_opt_;     /**< optimality gap tolerance */
	double eps_e_;       /**< linearization error tolerance */
	double eps_g_;       /**< gradient tolerance */
	double tau_min_;     /**< minimum bundle penalty parameter */
	double tau_;         /**< bundle penalty parameter */
	double upperbound_;  /**< upper bound */
	double valueAtProx_; /**< objective value at proximal point */
	double level_;       /**< level parameter */
	double gg_;          /**< 2-norm of subgradient */

	bool isSolved_; /**< Is the problem solved? */

	int nrows_;     /**< number of original rows (excluding row cuts) */
	int ncols_;     /**< number of original columns (excluding column cuts) */
	int nzcnt_;     /**< number of nonzero elements (excluding cuts) */

	double modelObjval_; /**< objective function value without proximal term */
	vector<bool> activeCols_;         /**< indication of active columns */
	vector<double*> lhsCouplingRows_; /**< LHS of coupling rows from previous iterations */
	vector<double>  dynObjs_;         /**< objective function coefficient dynamically generated */
	vector<int> dynids_;              /**< subproblem ids */

	/** pre-calculation for Hessian */
	double bb_;   /**< squared rhs of the coupling constraints */
};

#endif /* SRC_SOLVER_DUALDECOMP_DDMASTERDSB_H_ */
