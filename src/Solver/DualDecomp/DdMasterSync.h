/*
 * DdMasterSync.h
 *
 *  Created on: Mar 25, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DUALDECOMP_DDMASTERSYNC_H_
#define SRC_SOLVER_DUALDECOMP_DDMASTERSYNC_H_

#include <Solver/DualDecomp/DdMaster.h>

class DdMasterSync: public DdMaster {

	friend class DdMWSync;

public:

	/** constructor */
	DdMasterSync(
			DspParams *  par,     /**< parameter pointer */
			DecModel *   model,   /**< model pointer */
			DspMessage * message, /**< message pointer */
			int nworkers          /**< number of workers */);

	/** destructor */
	virtual ~DdMasterSync();

	/** initialize */
	virtual DSP_RTN_CODE init();

	/** update problem */
	virtual DSP_RTN_CODE updateProblem() {return DSP_RTN_OK;}

protected:

	int worker_;            /**< worker ID in communication */
	int nsubprobs_;         /**< number of subproblems for the current worker */
	int * subindex_;        /**< array of subproblem indices */
	double * subprimobj_;   /**< subproblem primal objective values */
	double * subdualobj_;   /**< subproblem dual objective values */
	double ** subsolution_; /**< subproblem solution */
};

#endif /* SRC_SOLVER_DUALDECOMP_DDMASTERSYNC_H_ */
