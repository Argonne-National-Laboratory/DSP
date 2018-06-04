/*
 * DspBranchObj.h
 *
 *  Created on: Oct 11, 2016
 *      Author: Kibaek Kim
 */

#ifndef SRC_TREESEARCH_DSPBRANCHOBJ_H_
#define SRC_TREESEARCH_DSPBRANCHOBJ_H_

#include <vector>
#include <CoinPackedVector.hpp>

class DspBranchObj {
private:
	//std::vector<int> index_;
	std::vector<CoinPackedVector*> vecs_; /** disjunction vector */
	std::vector<double> lb_;
	std::vector<double> ub_;

public:
	double bestBound_; /**< best bound */
	std::vector<double> dualsol_;
	int direction_;
	/** The solution estimate. The smaller the better. */
	double solEstimate_;

public:
	DspBranchObj(): bestBound_(-COIN_DBL_MAX), direction_(1), solEstimate_(0.0) {}
	virtual ~DspBranchObj() {
		for (auto it = vecs_.begin(); it != vecs_.end(); it++)
			delete *it;
	}

	void push_back(CoinPackedVector* vec, double lb, double ub) {
		vecs_.push_back(vec);
		lb_.push_back(lb);
		ub_.push_back(ub);
		vec = NULL;
	}

	void push_back(int index, double lb, double ub) {
		double coef = 1.0;
		push_back(new CoinPackedVector(1, &index, &coef), lb, ub);
	}

	/** number of branching objects */
	int getNumObjs() const {return vecs_.size();}

	/** get branching index for variable branching */
	int getIndex(int j) const {return vecs_[j]->getIndices()[0];}

	/** get branching vector for disjunction branching */
	const CoinPackedVector* getVector(int j) const {return vecs_[j];}

	double getLb(int j) const {return lb_[j];}
	double getUb(int j) const {return ub_[j];}
};

#endif /* SRC_TREESEARCH_DSPBRANCHOBJ_H_ */
