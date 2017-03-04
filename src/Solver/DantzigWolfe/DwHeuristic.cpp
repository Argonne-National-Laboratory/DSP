/*
 * DwHeuristic.cpp
 *
 *  Created on: Feb 12, 2017
 *      Author: kibaekkim
 */

#include <memory>
#include <DantzigWolfe/DwHeuristic.h>
#include <DantzigWolfe/DwMaster.h>
#include <DantzigWolfe/DwModel.h>

DwRounding::DwRounding(std::string name, DspModel &model):
DspHeuristic(name, model) {}

int DwRounding::solution(double &objective, std::vector<double> &solution) {

	/** create branching objects */
	std::shared_ptr<DspBranch> branch(new DspBranch);

	int found = 0;
	double stime = CoinGetTimeOfDay();
	std::shared_ptr<DwMaster> master(dynamic_cast<DwMaster*>(model_->getSolver()->clone()));
	//printf("Time to copy DwMaster: %.10f seconds\n", CoinGetTimeOfDay() - stime);
	std::vector<double> primsol(model_->getPrimalSolution());

	/** round and fix */
	double rounded;
	for (int j = 0; j < master->ncols_orig_; ++j) {
		if (master->ctype_orig_[j] != 'C') {
			/** round */
			rounded = round(primsol[j]);
			rounded = std::min(rounded, master->cubd_node_[j]);
			rounded = std::max(rounded, master->clbd_node_[j]);
			/** fix */
			branch->index_.push_back(j);
			branch->lb_.push_back(rounded);
			branch->ub_.push_back(rounded);
		}
	}
	branch->bestBound_ = master->getBestDualObjective();
	branch->dualsol_.assign(master->getBestDualSolution(), master->getBestDualSolution() + master->nrows_);

	master->setBranchingObjects(branch.get());

	double gaptol = master->getParPtr()->getDblParam("DW/GAPTOL");
	master->getParPtr()->setDblParam("DW/GAPTOL", 0.0);

	DSP_RTN_CHECK_THROW(master->solve());

	master->getParPtr()->setDblParam("DW/GAPTOL", gaptol);

	switch (master->getStatus()) {
	case DSP_STAT_OPTIMAL:
	case DSP_STAT_FEASIBLE:
	case DSP_STAT_LIM_ITERorTIME: {
		if (master->getPrimalObjective() < objective) {
			objective = -master->getBestDualObjective();
			/** parse solution */
			int cpos = 0;
			solution.resize(master->ncols_orig_);
			std::fill(solution.begin(), solution.end(), 0.0);
			for (auto it = master->cols_generated_.begin(); it != master->cols_generated_.end(); it++) {
				if ((*it)->active_) {
					for (int i = 0; i < (*it)->x_.getNumElements(); ++i) {
						if ((*it)->x_.getIndices()[i] < master->ncols_orig_)
							solution[(*it)->x_.getIndices()[i]] += (*it)->x_.getElements()[i] * master->getPrimalSolution()[cpos];
					}
					cpos++;
				}
			}
			found = 1;
		}
		break;
	}
	default:
		break;
	}

	return found;
}
