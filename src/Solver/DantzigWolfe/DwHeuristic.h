/*
 * DwHeuristic.h
 *
 *  Created on: Feb 12, 2017
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DANTZIGWOLFE_DWHEURISTIC_H_
#define SRC_SOLVER_DANTZIGWOLFE_DWHEURISTIC_H_

#include "TreeSearch/DspHeuristic.h"

class DwRounding: public DspHeuristic {
public:
	DwRounding(std::string name, DspModel &model);
	~DwRounding() {}

	virtual int solution(double &objective, std::vector<double> &solution);
};

#endif /* SRC_SOLVER_DANTZIGWOLFE_DWHEURISTIC_H_ */
