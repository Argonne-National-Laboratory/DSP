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

/** check whether solution is duplicate or not */
inline bool duplicateVector(
		CoinPackedVector * vec,
		vector<CoinPackedVector*> vecs)
{
	bool dup = false;

	/** number of saved solutions */
	int num = vecs.size();
	DSPdebugMessage2("number of vectors %d\n", num);
	for (int i = num - 1; i >= 0; --i)
	{
#ifdef DSP_DEBUG2
		DSPdebugMessage("vecs[%d]:\n", i);
		DspMessage::printArray(vecs[i]);
#endif
		if (vec->getNumElements() != vecs[i]->getNumElements() ||
			vec->getMinIndex() != vecs[i]->getMinIndex() ||
			vec->getMaxIndex() != vecs[i]->getMaxIndex() ||
			fabs(vec->infNorm() - vecs[i]->infNorm()) > 1.0e-8 ||
			fabs(vec->oneNorm() - vecs[i]->oneNorm()) > 1.0e-8 ||
			fabs(vec->sum() - vecs[i]->sum()) > 1.0e-8 ||
			fabs(vec->twoNorm() - vecs[i]->twoNorm()) > 1.0e-8)
			continue;
		if (vec->isEquivalent(*vecs[i]))
		{
			dup = true;
			break;
		}
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
