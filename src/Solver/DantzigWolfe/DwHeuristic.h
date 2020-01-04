/*
 * DwHeuristic.h
 *
 *  Created on: Feb 12, 2017
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DANTZIGWOLFE_DWHEURISTIC_H_
#define SRC_SOLVER_DANTZIGWOLFE_DWHEURISTIC_H_

#include "TreeSearch/DspHeuristic.h"
#include "Model/TssModel.h"
#include "Solver/DantzigWolfe/DwMaster.h"

class DwRounding: public DspHeuristic {
public:
	DwRounding(std::string name, DspModel &model):
		DspHeuristic(name, model) {}
	~DwRounding() {}

	virtual int solution(double &objective, std::vector<double> &solution);
};

class DwSmip: public DspHeuristic {
public:
	DwSmip(std::string name, DspModel &model):
		DspHeuristic(name, model) {}
	~DwSmip() {}

	virtual int solution(double &objective, std::vector<double> &solution);
	virtual void fixSolution(TssModel* tss, DwMaster* master, DspBranchObj* branch, std::vector<double>& sol);
};

#endif /* SRC_SOLVER_DANTZIGWOLFE_DWHEURISTIC_H_ */
