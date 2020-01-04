/*
 * DspHeuristic.h
 *
 *  Created on: Feb 13, 2017
 *      Author: kibaekkim
 */

#ifndef SRC_TREESEARCH_DSPHEURISTIC_H_
#define SRC_TREESEARCH_DSPHEURISTIC_H_

#include <vector>
#include <string>

class DspModel;

class DspHeuristic {
public:
	DspHeuristic(std::string name, DspModel &model) {
		name_ = name;
		model_ = &model;
	}

	virtual ~DspHeuristic() {}

	/** returns 0 if no solution, 1 if valid solution
	 * with better objective value than one passed in */
	virtual int solution(double &objective, std::vector<double> &solution) = 0;

	const char* name() {return name_.c_str();}

protected:
	std::string name_;
	DspModel *model_;
};

#endif /* SRC_TREESEARCH_DSPHEURISTIC_H_ */
