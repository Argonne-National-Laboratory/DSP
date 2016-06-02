/*
 * DdMasterAtr.h
 *
 *  Created on: Mar 22, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DUALDECOMP_DDMASTERATR_H_
#define SRC_SOLVER_DUALDECOMP_DDMASTERATR_H_

#include <Solver/DualDecomp/DdMasterTr.h>

class DdMasterAtr: public DdMasterTr {

	friend class DdMWAsync;

public:

	/** constructor */
	DdMasterAtr(
			DspParams *  par,
			DecModel *   model,
			DspMessage * message,
			int          nworkers,
			int          maxnumsubprobs);

	/** desctructor */
	virtual ~DdMasterAtr();

	/** initialize */
	virtual DSP_RTN_CODE init();

	/** solve */
	virtual DSP_RTN_CODE solve();

	/** update problem */
	virtual DSP_RTN_CODE updateProblem();

	/** termination test */
	DSP_RTN_CODE terminationTest();

private:

	vector<int> worker_;
	vector<int> nsubprobs_;    /**< number of subproblems for the current worker */
	vector<int*> subindex_;        /**< array of subproblem indices */
	vector<double*> subprimobj_;   /**< subproblem primal objective values */
	vector<double*> subdualobj_;   /**< subproblem dual objective values */
	vector<double**> subsolution_; /**< subproblem solution */

private:

	int *     nlastcuts_;         /**< number of cuts generated at the last iteration */
	double ** primsol_to_worker_; /**< primal solution (theta and lambda) given to each worker */
};

#endif /* SRC_SOLVER_DUALDECOMP_DDMASTERATR_H_ */
