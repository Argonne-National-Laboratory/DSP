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

public:
	/** default constructor */
	DspNodeSolution():
		AlpsSolution(),
		objvalue_(COIN_DBL_MAX) {}

	DspNodeSolution(std::vector<double>& sol, double obj):
		AlpsSolution(),
		solution_(sol),
		objvalue_(obj) {}

	/** default destructor */
	virtual ~DspNodeSolution() {}

	virtual void print(std::ostream& os) const {
		os << "Objective value = " << objvalue_ << std::endl;
		for (unsigned j = 0; j < solution_.size(); ++j) {
			if (fabs(solution_[j]) > 1.0e-8)
				os << "x[" << j << "] = " << solution_[j] << std::endl;
		}
	}

	std::vector<double> solution_; /**< solution */
	double objvalue_;  /**< objective value */
};

#endif /* SRC_TREESEARCH_DSPNODESOLUTION_H_ */
