/*
 * DspUtility.h
 *
 *  Created on: Jul 13, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_UTILITY_DSPUTILITY_H_
#define SRC_UTILITY_DSPUTILITY_H_

#include <vector>
/** Coin */
#include "OsiSolverInterface.hpp"
/** Dsp */
#include "Utility/DspRtnCodes.h"
#include "Utility/DspMessage.h"

using namespace std;

inline bool myduplicatetolerance(double i, double j) {
	return (fabs(i-j) < 1.0e-8);
}

/** check whether solution is duplicate or not */
inline bool duplicateVector(
		CoinPackedVector * vec,
		vector<CoinPackedVector*> vecs)
{
	bool dup = false;

	// printf("vec (%d):\n", vec->getNumElements());
	// DspMessage::printArray(vec);
	/** number of saved solutions */
	int num = vecs.size();
	// printf("number of vectors %d\n", num);
	
	for (int i = num - 1; i >= 0; --i)
	{
		// printf("vecs[%d] (%d):\n", i, vecs[i]->getNumElements());
		// DspMessage::printArray(vecs[i]);
		// printf("%sequal size\n", vec->getNumElements() == vecs[i]->getNumElements() ? "" : "not ");
		// printf("%d size\n", vec->getNumElements());
		// printf("%sequal indices\n", std::equal(vec->getIndices(), vec->getIndices()+vec->getNumElements(), vecs[i]->getIndices()) ? "" : "not ");
		// printf("%sequal elements\n", std::equal(vec->getElements(), vec->getElements()+vec->getNumElements(), vecs[i]->getElements(), myduplicatetolerance) ? "" : "not ");
		// printf("%s\n", dup ? "duplicate" : "new");

		if (vec->getNumElements() == vecs[i]->getNumElements()) {
			if (vec->getNumElements() == 0)
				dup = true;
			else if (std::equal(vec->getIndices(), vec->getIndices()+vec->getNumElements(), vecs[i]->getIndices()) && 
				std::equal(vec->getElements(), vec->getElements()+vec->getNumElements(), vecs[i]->getElements(), myduplicatetolerance))
				dup = true;
		}
		// printf("%s\n", dup ? "duplicate" : "new");
		if (dup) break;
	}

	return dup;
}

/** convert coin-status to dsp-status */
inline void convertCoinToDspStatus(const OsiSolverInterface* si, int& status) {
	if (si->isProvenOptimal())
		status = DSP_STAT_OPTIMAL;
	else if (si->isProvenPrimalInfeasible())
		status = DSP_STAT_PRIM_INFEASIBLE;
	else if (si->isProvenDualInfeasible())
		status = DSP_STAT_DUAL_INFEASIBLE;
	else if (si->isPrimalObjectiveLimitReached())
		status = DSP_STAT_LIM_PRIM_OBJ;
	else if (si->isDualObjectiveLimitReached())
		status = DSP_STAT_LIM_DUAL_OBJ;
	else if (si->isIterationLimitReached())
		status = DSP_STAT_LIM_ITERorTIME;
	else if (si->isAbandoned())
		status = DSP_STAT_ABORT;
	else
		status = DSP_STAT_UNKNOWN;
}


#endif /* SRC_UTILITY_DSPUTILITY_H_ */
