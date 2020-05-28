/*
 * DdMasterTr.h
 *
 *  Created on: Feb 12, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DUALDECOMP_DDMASTERTR_H_
#define SRC_SOLVER_DUALDECOMP_DDMASTERTR_H_

#include "Solver/DualDecomp/DdMaster.h"

/** A class for implementing the trust-region master solver. */
class DdMasterTr: public DdMaster {
public:

	/** A default constructor. */
	DdMasterTr(
			DecModel *   model,   /**< model pointer */
			DspParams *  par,     /**< parameter pointer */
			DspMessage * message /**< message pointer */);

	/** A copy constructor. */
	DdMasterTr(const DdMasterTr& rhs);

	/** A default desctructor. */
	virtual ~DdMasterTr();

	/** A clone function */
	virtual DdMasterTr* clone() const {
		return new DdMasterTr(*this);
	}

	/** A virtual member for initializing solver. */
	virtual DSP_RTN_CODE init();

	/** A virtual member for solving problem. */
	virtual DSP_RTN_CODE solve();

	/** update problem */
	virtual DSP_RTN_CODE updateProblem();

	/** termination test */
	virtual DSP_RTN_CODE terminationTest();

protected:

	/** create problem */
	virtual DSP_RTN_CODE createProblem();

	/** is solution trust region boundary? */
	virtual bool isSolutionBoundary(double eps = 1.0e-6);

	/** add cuts */
	virtual int addCuts(
			bool possiblyDel = true /**< possibly delete cuts*/);

	/** possibly delete cuts */
	virtual DSP_RTN_CODE possiblyDeleteCuts(
			double subobjval /**< sum of subproblem objective values */);

	/** possibly delete cuts */
	virtual DSP_RTN_CODE possiblyDeleteCutsOsi(
			double subobjval /**< sum of subproblem objective values */);

	/** recruite cuts */
	virtual int recruiteCuts();

	/** remove all cuts */
	virtual DSP_RTN_CODE removeAllCuts();

	/** change trust region */
	virtual DSP_RTN_CODE setTrustRegion(double stability_param, double * stability_center);

protected:

	int nthetas_;  /**< number of thetas */
	int nlambdas_; /**< number of lambdas */
	int nus_;      /**< DRO: number of u's */
	int nPs_;      /**< DRO: number of Ps */

	double ** obj_reco_; /**< DRO: recourse objective coefficient for each scenario */

	double stability_param_;    /**< stability parameter */
	double * stability_center_; /**< stability center */
	int trcnt_;                 /**< trust region counter */
	int numIters_;              /**< number of iterations */
	double cputime_elapsed_;    /**< cpu time elapsed */
	double walltime_elapsed_;   /**< wall time elapsed */

	bool isSolved_; /**< indicating whether problem is ever solved */


	/** TODO: The following could be encapsulated by a class. */
	OsiCuts *      cuts_;           /**< cut pool (can have column cuts and row cuts) */
	vector<int>    cuts_age_;       /**< ages of cuts: -1 (not added) */
	vector<bool>   possiblyDelete_; /**< indicating whether cuts are possibly deleted. */
	vector<double> masterobjsAtCutAdd_; /**< master objective values when cuts were generated */
	int            ncuts_minor_;    /**< number of cuts generated at minor iterations */
	double         cutdel_param_;   /**< cut deletion parameter */
	double         linerr_;         /**< linearization error */

	/** parameters */
	bool parTr_;            /**< enable/disable trust region */
	double parTrSize_;      /**< trust region size */
	bool parTrDecrease_;    /**< enable decreasing trust region */
	int parNumCutsPerIter_; /**< number of cuts added per iteration (determines dimension of theta) */
	int parMasterAlgo_;     /**< algorithm for solving master problem */
	int parLogLevel_;       /**< display level */
	int nstalls_;           /**< number of stalling iterations */
};

#endif /* SRC_SOLVER_DUALDECOMP_DDMASTERTR_H_ */
