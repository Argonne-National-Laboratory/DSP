/*
 * DspUtility.h
 *
 *  Created on: Jul 13, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_UTILITY_DSPUTILITY_H_
#define SRC_UTILITY_DSPUTILITY_H_

#include <vector>
#include "CoinHelperFunctions.hpp"
#include "Utility/DspMessage.h"

using namespace std;

/** check whether solution is duplicate or not */
bool duplicateVector(
		CoinPackedVector * vec,
		vector<CoinPackedVector*> vecs)
{
	bool dup = false;

	/** number of saved solutions */
	int num = vecs.size();
	DSPdebugMessage("number of vectors %d\n", num);
	for (int i = num - 1; i >= 0; --i)
	{
#ifdef DSP_DEBUG
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


#endif /* SRC_UTILITY_DSPUTILITY_H_ */
