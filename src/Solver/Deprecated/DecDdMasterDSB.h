/*
 * DecDdMasterDSB.h
 *
 *  Created on: Aug 12, 2015
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DECDDMASTERDSB_H_
#define SRC_SOLVER_DECDDMASTERDSB_H_

#include <Solver/DecDdMaster.h>
#include "SolverInterface/SolverInterfaceOoqp.h"

using namespace std;

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
class DecDdMasterDSB: public DecDdMaster {
public:

	/** default constructor */
	DecDdMasterDSB(DspParams * par) :
		DecDdMaster(par),
		si_(NULL),
		model_(NULL),
		prox_(NULL), lambda_(NULL),
		phi_t_(10.0), phi_l_(1000.0),
		alpha_t_(0.1), alpha_l_(0.1),
		eps_opt_(1.0e-6), eps_e_(1.0e-6), eps_g_(1.0e-6),
		tau_min_(1.0e-6), tau_(1.0),
		upperbound_(1.0e+20), valueAtProx_(-1.0e+20), level_(0.0), gg_(0.0),
		isSolved_(false),
		nrows_(0), ncols_(0), nzcnt_(0), ncoupling_(0), nsubprobs_(0), ncols0_(0),
		modelObjval_(0.0),
		bb_(0.0)
	{
		/** nothing to do */
	}

	/** default destructor */
	virtual ~DecDdMasterDSB();

	/** create problem */
	virtual DSP_RTN_CODE createProblem(DecModel * model);

	/** update problem: may update dual bound */
	virtual DSP_RTN_CODE updateProblem(
			double primal_bound,     /**< primal bound of the original problem */
			double & dual_bound,     /**< dual bound of the original problem */
			double * primal_objvals, /**< objective values of subproblems */
			double * dual_objvals,   /**< objective values of subproblems */
			double ** solution       /**< subproblem solutions */);

	/** solve problem */
	virtual DSP_RTN_CODE solve();

	/** solution status */
	virtual DSP_RTN_CODE getStatus() {return si_->getStatus();}

	/** get Lagrangian multiplier */
	virtual const double * getLagrangian() {return lambda_;}

	/** termination test */
	virtual bool terminate();

private:

	/** manage bundle */
	DSP_RTN_CODE manageBundle(
			double * objvals /**< subproblem objective values */);

	/** apply level change */
	DSP_RTN_CODE applyLevelChange();

	/** update dynamic objective coefficients */
	DSP_RTN_CODE updateDynObjs();

public:

	SolverInterfaceOoqp * si_; /**< solver interface */

private:

	DecModel * model_; /**< model */

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
	int ncoupling_; /**< number of coupling constraints */
	int nsubprobs_; /**< number of scenario subproblems */
	int ncols0_;    /**< number of fist-stage columns (stochastic only) */

	double modelObjval_; /**< objective function value without proximal term */
	vector<bool> activeCols_;         /**< indication of active columns */
	vector<double*> lhsCouplingRows_; /**< LHS of coupling rows from previous iterations */
	vector<double>  dynObjs_;         /**< objective function coefficient dynamically generated */
	vector<int> dynids_;              /**< subproblem ids */

	/** pre-calculation for Hessian */
	double bb_;   /**< squared rhs of the coupling constraints */
};

#endif /* SRC_SOLVER_DECDDMASTERDSB_H_ */
