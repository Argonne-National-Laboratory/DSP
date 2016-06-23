/*
 * DdWorkerCGBd.h
 *
 *  Created on: May 23, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DUALDECOMP_DDWORKERCGBD_H_
#define SRC_SOLVER_DUALDECOMP_DDWORKERCGBD_H_

#include "Solver/DualDecomp/DdWorkerCG.h"
#include "Solver/Benders/BdSub.h"

class DdWorkerCGBd: public DdWorkerCG {

public:

	enum {
		Opt = 0,
		Feas,
		None
	};

	/** constructor */
	DdWorkerCGBd(DspParams * par, DecModel * model, DspMessage * message);

	/** destructor */
	virtual ~DdWorkerCGBd();

	/** initialize */
	virtual DSP_RTN_CODE init();

	/** generate cuts */
	virtual DSP_RTN_CODE generateCuts(
			CoinPackedVector * solutions, /**< solutions to evaluate */
			OsiCuts *   cuts,             /**< cuts */
			vector<int> & cuttype         /**< cut type: Opt or Feas */);

	/** get type of worker */
	virtual int getType() {return CGBd;}

protected:

	BdSub * bdsub_; /**< pointer to Benders subproblem */
};

#endif /* SRC_SOLVER_DUALDECOMP_DDWORKERCGBD_H_ */
