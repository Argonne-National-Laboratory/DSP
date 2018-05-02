/*
 * DspBranchObj.h
 *
 *  Created on: Oct 11, 2016
 *      Author: Kibaek Kim
 */

#ifndef SRC_TREESEARCH_DSPBRANCHOBJ_H_
#define SRC_TREESEARCH_DSPBRANCHOBJ_H_

#include <vector>

struct DspBranchObj {
	std::vector<int> index_;
	std::vector<double> lb_;
	std::vector<double> ub_;
	double bestBound_; /**< best bound */
	std::vector<double> dualsol_;
	int direction_;

	void push_back(int index, double lb, double ub) {
		index_.push_back(index);
		lb_.push_back(lb);
		ub_.push_back(ub);
	}
};

#endif /* SRC_TREESEARCH_DSPBRANCHOBJ_H_ */
