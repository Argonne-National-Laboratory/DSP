/*
 * DspHeuristic.h
 *
 *  Created on: Feb 13, 2017
 *      Author: kibaekkim
 */

#ifndef SRC_TREESEARCH_DSPHEURISTIC_H_
#define SRC_TREESEARCH_DSPHEURISTIC_H_

#include "Solver/DecSolver.h"

class DspHeuristic {
public:
	DspHeuristic(std::string name, DecSolver &solver) {
		name_ = name;
		solver_ = &solver;
	}

	virtual ~DspHeuristic() {}

	/** returns 0 if no solution, 1 if valid solution
	 * with better objective value than one passed in */
	virtual int solution(double &objective, std::vector<double> &solution) = 0;

	const char* name() {return name_.c_str();}

protected:
	std::string name_;
	DecSolver *solver_;
};

#endif /* SRC_TREESEARCH_DSPHEURISTIC_H_ */
