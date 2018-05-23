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
#include "Model/TssModel.h"

int DwRounding::solution(double &objective, std::vector<double> &solution) {

	/** create branching objects */
	std::shared_ptr<DspBranchObj> branch(new DspBranchObj);

	int found = 0;
	double stime = CoinGetTimeOfDay();
	std::shared_ptr<DwMaster> master(dynamic_cast<DwMaster*>(model_->getSolver()->clone()));
	//printf("Time to copy DwMaster: %.10f seconds\n", CoinGetTimeOfDay() - stime);

	/** round and fix */
	double rounded;
	for (int j = 0; j < master->ncols_orig_; ++j) {
		if (master->ctype_orig_[j] != 'C') {
			/** round */
			rounded = round(model_->getPrimalSolution()[j]);
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
			objective = master->getBestDualObjective();
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

int DwSmip::solution(double &objective, std::vector<double> &solution) {

	int found = 0;
	TssModel* tss = NULL;
	DecSolver* solver = model_->getSolver();
	DecModel* model = solver->getModelPtr();
	DspMessage* message = solver->getMessagePtr();

	/** get stochastic model pointer */
	if (model->isStochastic()) {
		tss = dynamic_cast<TssModel*>(model);
	} else {
		return found;
	}

	/** create branching objects */
	std::shared_ptr<DspBranchObj> branch(new DspBranchObj);

	double stime = CoinGetTimeOfDay();
	std::shared_ptr<DwMaster> master(dynamic_cast<DwMaster*>(solver->clone()));
	//printf("Time to copy DwMaster: %.10f seconds\n", CoinGetTimeOfDay() - stime);

	/** round and fix */
	double fixed;
	for (int j = 0; j < tss->getNumCols(0); ++j) {
		fixed = model_->getPrimalSolution()[j];
		if (master->ctype_orig_[j] != 'C') {
			/** round */
			double rounded = round(fixed);
			fixed = std::max(std::min(rounded, master->cubd_node_[j]), master->clbd_node_[j]);
			message->print(5, "fixing variable %d to %e\n", j, fixed);
		}
		/** fix */
		for (int s = 0; s < tss->getNumScenarios(); ++s)
			branch->push_back(tss->getNumCols(0) * s + j, fixed, fixed);
	}
	//branch->bestBound_ = master->getBestDualObjective();
	//branch->dualsol_.assign(master->getBestDualSolution(), master->getBestDualSolution() + master->nrows_);

	master->setBranchingObjects(branch.get());

	int loglevel = master->getParPtr()->getIntParam("LOG_LEVEL");
	double gaptol = master->getParPtr()->getDblParam("DW/GAPTOL");
	master->getParPtr()->setIntParam("LOG_LEVEL", 0);
	master->getParPtr()->setDblParam("DW/GAPTOL", 0.0);

	DSP_RTN_CHECK_THROW(master->solve());

	master->getParPtr()->setIntParam("LOG_LEVEL", loglevel);
	master->getParPtr()->setDblParam("DW/GAPTOL", gaptol);

	switch (master->getStatus()) {
	case DSP_STAT_OPTIMAL:
	case DSP_STAT_FEASIBLE:
	case DSP_STAT_LIM_ITERorTIME: {
		if (master->getDualObjective() < objective) {
			objective = master->getBestDualObjective();
			message->print(1, "found a better upper bound %e.\n", objective);
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

