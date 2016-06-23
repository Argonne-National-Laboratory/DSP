/*
 * DecDdPrimalMaster.h
 *
 *  Created on: Dec 11, 2014
 *      Author: kibaekkim, ctjandra
 */

#ifndef SRC_SOLVER_DECDDPRIMALMASTER_H_
#define SRC_SOLVER_DECDDPRIMALMASTER_H_

#include "Solver/DecDdMaster.h"

class DecDdPrimalMaster: public DecDdMaster
{
public:

	/** constructor */
	DecDdPrimalMaster(DspParams * par) :
		DecDdMaster(par),
		si_(NULL),
		model_(NULL),
		nrows_(0),
		ncols_(0),
		ncoupling_(0),
		nsubprobs_(0),
		nCutsPerIter_(par->getIntParam("DD/NUM_CUTS_PER_ITER")),
		cuts_(NULL),
		ncuts_minor_(0),
		cutdel_param_(0.5),
		prox_(NULL),
		rho_(par->getDblParam("DD/TR/SIZE")),
		trcnt_(0),
		ubCuts_(NULL),
		isSolved_(false)
	{
		/** nothing to do */
	}

	/** default destructor */
	virtual ~DecDdPrimalMaster();

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

	virtual DSP_RTN_CODE getStatus() {return si_->getStatus();}

	/** get Lagrangian multiplier */
	virtual const double * getLagrangian() {return (solution_ + nCutsPerIter_);}

	/** termination test */
	virtual bool terminate();

	/** is solution trust region boundary? */
	virtual bool isSolutionBoundary(double eps = 1.0e-5);

	/** get rho */
	virtual double getRho() const {return rho_;}

	/** get proximal point */
	virtual const double * getProxPoint() const {return prox_;}

	/** set rho */
	virtual void setRho(double rho) {rho_ = rho;}

	/** change trust region */
	virtual void setTrustRegion(double rho, double * prox = NULL)
	{
		if (prox == NULL)
			prox = prox_;
		for (int j = nCutsPerIter_; j < ncols_; ++j)
		{
			double clbd = prox[j - nCutsPerIter_] - rho;
			double cubd = prox[j - nCutsPerIter_] + rho;
			if (model_->getSenseCouplingRow(j - nCutsPerIter_) == 'L')
				clbd = CoinMax((double) 0.0, clbd); /* lambda >= 0 */
			else if (model_->getSenseCouplingRow(j - nCutsPerIter_) == 'G')
				cubd = CoinMin((double) 0.0, cubd); /* lambda <= 0 */
			si_->setColBounds(j, clbd, cubd);
		}
	}

	/** set print level */
	virtual void setPrintLevel(int level) {si_->setPrintLevel(level);}

private:

	/** aging cuts */
	virtual DSP_RTN_CODE agingCuts();

	/** add cuts */
	virtual int addCuts(
			double * objvals,            /**< objective values of subproblems */
			double ** solution,          /**< subproblem solutions */
			bool      possiblyDel = true /**< possibly delete cuts*/);

	/** possibly delete cuts */
	virtual DSP_RTN_CODE possiblyDeleteCuts(
			double sub_objval /**< subproblem objective value */);

	/** recruit cuts: investigating and adding effective cuts (that were deleted) in the cut pool */
	virtual int recruiteCuts();

	/** remove all cuts from solver interface */
	virtual DSP_RTN_CODE removeAllCuts();

public:

	SolverInterface * si_; /**< solver interface */

private:

	DecModel * model_; /**< model */

	int nrows_;       /**< number of original rows (excluding row cuts) */
	int ncols_;       /**< number of original columns (excluding column cuts) */
	int ncoupling_;   /**< number of coupling constraints */
	int nsubprobs_;   /**< number of scenario subproblems */
	int nCutsPerIter_;/**< number of cuts added per iteration */

	OsiCuts *      cuts_;           /**< cut pool (can have column cuts and row cuts) */
	vector<int>    cuts_age_;       /**< ages of cuts: -1 (not added) */
	vector<bool>   possiblyDelete_; /**< indicating whether cuts are possibly deleted. */
	vector<double> master_objvals_; /**< master objective values when cuts were generated */
	int            ncuts_minor_;    /**< number of cuts generated at minor iterations */
	double         cutdel_param_;   /**< cut deletion parameter */

	double * prox_;       /**< proximal point for trust region */
	double rho_;          /**< trust region size */
	int trcnt_;           /**< trust region counter */

	OsiCuts * ubCuts_; /**< upper bound cuts for regularization */

	bool isSolved_;   /**< indicating whether problem is ever solved */
};

#endif /* SRC_SOLVER_DECDDPRIMALMASTER_H_ */
