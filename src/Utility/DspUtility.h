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
	DSPdebugMessage2("number of vectors %d\n", num);
	for (int i = num - 1; i >= 0; --i)
	{
#ifdef DSP_DEBUG2
		DSPdebugMessage("vecs[%d]:\n", i);
		DspMessage::printArray(vecs[i]);
#endif
		if (vec->getNumElements() == 0 || vecs[i]->getNumElements() == 0) {
		  if (vec->getNumElements() == 0 && vecs[i]->getNumElements() == 0)
		    dup = true;
		  else
		    dup = false;
		  break;
		} else {
		  dup = (vec->getNumElements() == vecs[i]->getNumElements() &&
		    std::equal(vec->getIndices(),vec->getIndices()+vec->getNumElements(),
			     vecs[i]->getIndices()) &&
		    std::equal(vec->getElements(),vec->getElements()+vec->getNumElements(),
			       vecs[i]->getElements()));
		}
	}

	return dup;
}


#endif /* SRC_UTILITY_DSPUTILITY_H_ */
