/*
 * BdSub.h
 *
 *  Created on: Feb 18, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_BENDERS_BDSUB_H_
#define SRC_SOLVER_BENDERS_BDSUB_H_

#include <assert.h>
#include "Utility/DspMacros.h"
#include "Utility/DspRtnCodes.h"
#include "Utility/DspParams.h"
#include "SolverInterface/DspOsi.h"
#include "Model/DecModel.h"

/**
 * This class defines a Benders cut generator.
 */
class BdSub {
public:

	/** A default constructor. */
	BdSub(DspParams * par);

	/** A copy constructor. */
	BdSub(const BdSub& rhs);

	/** A default destructor. */
	virtual ~BdSub();

	/** A clone function. */
	virtual BdSub* clone() const {
		return new BdSub(*this);
	}

	/** set subproblem indices */
	DSP_RTN_CODE setSubIndices(int size, const int * indices);

	/** load problem */
	DSP_RTN_CODE loadProblem(DecModel * model);

	/** generate Benders cut in raw format (without constructing) */
	int generateCuts(
			int            ncols,  /**< number of first-stage variables */
			const double * x,      /**< first-stage solution */
			double **      cutval, /** dense cut coefficients for each subproblem */
			double *       cutrhs  /** cut rhs for each subproblem */);

	/** evaluate recourse function */
	DSP_RTN_CODE evaluateRecourse(
		const double * x, /**< [in] first-stage solution */
		double * objvals  /**< [out] objective values */
	);

public:

	/** get objective value */
	double getObjective(int i) {return objvals_[i];}

	/** get solution */
	const double * getSolution(int i) {return solutions_[i];}

	/** get status for subproblem subindices_[i] */
	DSP_RTN_CODE getStatus(int i) {return status_[i];}

	/** get statuses */
	const int* getStatuses() {return status_;}

	/** get number of subproblems */
	int getNumSubprobs() {return nsubprobs_;}

	/** get subproblem index */
	int getSubprobIndex(int i) {return subindices_[i];}

	/** get subproblem indices */
	int * getSubprobIndices() {return subindices_;}

	/** get number of columns */
	int getNumCols(int i) {return cglp_[i]->si_->getNumCols();}

	/** does recourse have integer variables? */
	bool has_integer() {return recourse_has_integer_;}

private:

	static DspOsi * createDspOsi(int solver);

	/** solve one subproblem. this is a body of loop in gutsOfGenerateCuts */
	static void solveOneSubproblem(
			BdSub *     cgl,
			int            s,                /**< scenario index */
			const double * x,                /**< first-stage solution */
			double **      Tx,               /**< Tx */
			double **      cutval,           /**< Benders cut body */
			double *       cutrhs,           /**< Benders cut RHS */
			int            enableOptCuts = 1 /**< whether to generate optimality cuts or not */);

	/** solve one integer subproblem */
	static DSP_RTN_CODE solveOneIntegerSubproblem(
			BdSub *     cgl,
			int            s,     /**< scenario index */
			const double * x,     /**< first-stage solution */
			double *       Tx,    /**< Tx */
			double *       objval /**< objective value */);

	/** solve feasibility problem */
	static DSP_RTN_CODE solveFeasProblem(
		OsiSolverInterface *si, /**< [in] subproblem solver interface */
		int &nAddedCols /**< [out] number of columns added */);

	/** change feasibility problem to original problem */
	static DSP_RTN_CODE chgToOrgProblem(
		OsiSolverInterface *si, /**< [in] subproblem solver interface */
		const double *obj,		/**< [in] original objective function */
		int &nAddedCols /**< [out] number of columns added */);

	/** calculate cut elements and rhs */
	static int calculateCutElements(
			int nrows,                         /**< [in] number of rows in subproblem */
			int ncols,                         /**< [in] number of columns in subproblem */
			const CoinPackedMatrix * mat_tech, /**< [in] technology matrix */
			const double * rlbd,               /**< [in] row lower bounds */
			const double * rubd,               /**< [in] row upper bounds */
			const double * clbd,               /**< [in] column lower bounds */
			const double * cubd,               /**< [in] column upper bounds */
			const double * pi,                 /**< [in] dual variables corresponding to constraints */
			const double * rc,                 /**< [in] reduced cost corresponding to variables */
			double *       cutval,             /**< [out] cut coefficients */
			double &       cutrhs              /**< [out] cut rhs */);

protected:

	DspParams * par_; /**< parameters */

	int   nsubprobs_;  /**< number of subproblems */
	int * subindices_; /**< subproblem indices indices for cut generation */
	double* probability_; /**< probability; 1.0 for non-stochastic model */

	CoinPackedMatrix   ** mat_mp_;     /**< array of matrix corresponding to master problem part */
	DspOsi ** cglp_;                   /**< array of Cut Generation LP */
	CoinWarmStart **      warm_start_; /**< warm start information for each subproblem */
	double *              objvals_;    /**< subproblem objective values */
	double **             solutions_;  /**< subproblem solutions */
	DSP_RTN_CODE *        status_;     /**< subproblem solution status */
	bool recourse_has_integer_;        /**< whether the recourse has integer variables */

	/** simple statistics */
	vector<string> names_statistics_;
	unordered_map<string, int> count_statistics_;
	unordered_map<string, double> time_statistics_;
};

#endif /* SRC_SOLVER_BENDERS_BDSUB_H_ */
