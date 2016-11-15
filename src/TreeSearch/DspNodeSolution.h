/*
 * DspNodeSolution.h
 *
 *  Created on: Oct 12, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_TREESEARCH_DSPNODESOLUTION_H_
#define SRC_TREESEARCH_DSPNODESOLUTION_H_

/** Coin */
#include "AlpsSolution.h"
/** Dsp */
#include "Utility/DspMacros.h"

class DspNodeSolution: public AlpsSolution {
protected:
	int len_;          /**< length of solution array */
	double* solution_; /**< solution */
	double objvalue_;  /**< objective value */

public:
	/** default constructor */
	DspNodeSolution():
		AlpsSolution(),
		len_(0),
		solution_(NULL),
		objvalue_(COIN_DBL_MAX) {}

	DspNodeSolution(int len, const double* sol, double obj):
		AlpsSolution(),
		len_(len),
		objvalue_(obj) {
		solution_ = new double [len_];
		CoinCopyN(sol, len, solution_);
	}

	/** default destructor */
	virtual ~DspNodeSolution() {
		FREE_ARRAY_PTR(solution_)
	}

	virtual void print(std::ostream& os) const {
	//	os << "Objective value = " << objvalue_ << std::endl;
		for (int j = 0; j < len_; ++j) {
			if (fabs(solution_[j]) > 1.0e-8)
				os << "x[" << j << "] = " << solution_[j] << std::endl;
		}
	}
};

#endif /* SRC_TREESEARCH_DSPNODESOLUTION_H_ */
