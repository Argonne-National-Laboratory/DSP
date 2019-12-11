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
#include "OsiCuts.hpp"

/** A Benders-cut generation class */
class DdWorkerCGBd: public DdWorkerCG {

public:

	enum {
		Opt = 0,
		Feas,
		None
	};

	/** A default constructor. */
	DdWorkerCGBd(
			DecModel *   model,  /**< model pointer */
			DspParams *  par,    /**< parameter pointer */
			DspMessage * message /**< message pointer */);

	/** A copy constructor. */
	DdWorkerCGBd(const DdWorkerCGBd& rhs);

	/** A default destructor. */
	virtual ~DdWorkerCGBd();

	/** A clone function */
	virtual DdWorkerCGBd* clone() const {
		return new DdWorkerCGBd(*this);
	}

	/** A virtual member to initialize the class. */
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
