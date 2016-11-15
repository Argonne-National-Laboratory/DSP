/*
 * DspNodeDesc.h
 *
 *  Created on: Sep 8, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_TREESEARCH_DSPNODEDESC_H_
#define SRC_TREESEARCH_DSPNODEDESC_H_

/** Coin */
#include "AlpsNodeDesc.h"
#include "CoinWarmStartBasis.hpp"
/** Dsp */
#include "Utility/DspMacros.h"
#include "TreeSearch/DspBranch.h"
#include "TreeSearch/DspModel.h"

/**
 * This class is a bridge to hand over the branching object to solver.
 */
class DspNodeDesc: public AlpsNodeDesc {
public:

	/** default constructor */
	DspNodeDesc() :
			AlpsNodeDesc(),
			branchdir_(1),
			branchobj_(NULL) {
		/** nothing to do */
	}

	/** constructor */
	DspNodeDesc(AlpsModel* m) :
			AlpsNodeDesc(m),
			branchdir_(1),
			branchobj_(NULL) {
		/** nothing to do */
	}

	/** constructor with branching object */
	DspNodeDesc(DspModel* m, int branchdir, DspBranch*& branchobj):
			AlpsNodeDesc(m),
			branchdir_(branchdir),
			branchobj_(branchobj) {
		branchobj = NULL;
	}

	/** default destructor */
	virtual ~DspNodeDesc() {
		FREE_PTR(branchobj_);
	}

	/** get branching object */
	const DspBranch* getBranchingObject() {return branchobj_;}

private:

	int branchdir_;        /**< branching direction: 1=up, -1=down */
	DspBranch* branchobj_; /**< branching object */
};

#endif /* SRC_TREESEARCH_DSPNODEDESC_H_ */
