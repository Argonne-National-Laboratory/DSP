/*
 * DwModel.cpp
 *
 *  Created on: Feb 13, 2017
 *      Author: kibaekkim
 */

#include <DantzigWolfe/DwModel.h>
#include <DantzigWolfe/DwHeuristic.h>

DwModel::DwModel(): DspModel() {}

DwModel::DwModel(DecSolver* solver): DspModel(solver) {
	heuristics_.push_back(new DwRounding("Rounding", *solver_));
}

DwModel::~DwModel() {
	for (unsigned i = 0; i < heuristics_.size(); ++i)
		delete heuristics_[i];
}

DSP_RTN_CODE DwModel::solve() {
	double bestprimobj = COIN_DBL_MAX;
	std::vector<double> solution;

	BGN_TRY_CATCH

	solver_->solve();

	switch (solver_->getStatus()) {
	case DSP_STAT_OPTIMAL:
	case DSP_STAT_FEASIBLE:
	case DSP_STAT_LIM_ITERorTIME:
		for (auto it = heuristics_.begin(); it != heuristics_.end(); it++) {
			printf("Running %s heuristic:\n", (*it)->name());
			int found = (*it)->solution(bestprimobj, solution);
			printf("found %d bestprimobj %+e\n", found, bestprimobj);
		}
		break;
	default:
		break;
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

