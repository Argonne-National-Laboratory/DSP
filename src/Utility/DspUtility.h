/*
 * DspUtility.h
 *
 *  Created on: Jul 13, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_UTILITY_DSPUTILITY_H_
#define SRC_UTILITY_DSPUTILITY_H_

//#define DSP_DEBUG
//#define DSP_DEBUG2

#include <vector>
#include "Utility/DspRtnCodes.h"
#include "Utility/DspMessage.h"

using namespace std;

inline bool is_integer(double v, double eps) {
	return fabs(v - floor(v + 0.5)) < eps;
}

inline bool myduplicatetolerance(double i, double j) {
	return (fabs(i-j) < 1.0e-8);
}

/** check whether solution is duplicate or not */
inline bool duplicateVector(
		CoinPackedVector * vec,
		vector<CoinPackedVector*> vecs)
{
	bool dup = false;

#ifdef DSP_DEBUG2
	DSPdebugMessage("vec (%d):\n", vec->getNumElements());
	DspMessage::printArray(vec);
#endif
	/** number of saved solutions */
	int num = vecs.size();
	// printf("number of vectors %d\n", num);
	
	for (int i = num - 1; i >= 0; --i)
	{
#ifdef DSP_DEBUG2
		DSPdebugMessage("vecs[%d] (%d):\n", i, vecs[i]->getNumElements());
		DspMessage::printArray(vecs[i]);
#endif
		if (vec->getNumElements() == vecs[i]->getNumElements()) {
			if (vec->getNumElements() == 0)
				dup = true;
			else if (std::equal(vec->getIndices(), vec->getIndices()+vec->getNumElements(), vecs[i]->getIndices()) && 
				std::equal(vec->getElements(), vec->getElements()+vec->getNumElements(), vecs[i]->getElements(), myduplicatetolerance))
				dup = true;
		}
		if (dup) break;
	}

	return dup;
}

#endif /* SRC_UTILITY_DSPUTILITY_H_ */
