/*
 * DwHeuristic.cpp
 *
 *  Created on: Feb 12, 2017
 *      Author: kibaekkim
 */

// #define DSP_DEBUG

#include <memory>
#include "Solver/DantzigWolfe/DwHeuristic.h"
#include "Solver/DantzigWolfe/DwModel.h"

int DwRounding::solution(double &objective, std::vector<double> &solution) {


	/** create branching objects */
	std::shared_ptr<DspBranchObj> branch(new DspBranchObj);

	int found = 0;
	std::shared_ptr<DwMaster> master(dynamic_cast<DwMaster*>(model_->getSolver()->clone()));
	//printf("Time to copy DwMaster: %.10f seconds\n", CoinGetTimeOfDay() - stime);
	DspMessage* message = master->getMessagePtr();

	/** primal solution */
	std::vector<double> refsol;
	dynamic_cast<DwModel*>(model_)->getRefSol(refsol);

	/** round and fix */
	double rounded;
	for (int j = 0; j < master->ncols_orig_; ++j) {
		if (master->ctype_orig_[j] != 'C') {
			/** round */
			rounded = round(refsol[j]);
			rounded = std::min(rounded, master->cubd_node_[j]);
			rounded = std::max(rounded, master->clbd_node_[j]);
			/** fix */
			branch->push_back(j, rounded, rounded);
		}
	}
	branch->bestBound_ = master->getBestDualObjective();
	branch->dualsol_.assign(master->getBestDualSolution(), master->getBestDualSolution() + master->nrows_);

	master->setBranchingObjects(branch.get());

	int loglevel = master->getParPtr()->getIntParam("LOG_LEVEL");
	int dwevalub = master->getParPtr()->getIntParam("DW/EVAL_UB");
	double gaptol = master->getParPtr()->getDblParam("DW/GAPTOL");
	message->logLevel_ = 0;
	master->getParPtr()->setIntParam("DW/EVAL_UB", -1);
	master->getParPtr()->setDblParam("DW/GAPTOL", 0.0);

	DSP_RTN_CHECK_THROW(master->solve());

	message->logLevel_ = loglevel;
	master->getParPtr()->setIntParam("DW/EVAL_UB", dwevalub);
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
#ifdef DSP_DEBUG2
	printf("objective = %e\n", objective);
	printf("solution (size %u):\n", solution.size());
	DspMessage::printArray(solution.size(), &solution[0]);
#endif

	int found = 0;
	TssModel* tss = NULL;
	DecSolver* solver = model_->getSolver();
	DwMaster* m = dynamic_cast<DwMaster*>(solver);
	DecModel* model = solver->getModelPtr();
	DspMessage* message = solver->getMessagePtr();

	/** get stochastic model pointer */
	if (model->isStochastic()) {
		tss = dynamic_cast<TssModel*>(model);
	} else {
		return found;
	}
	DSPdebugMessage("Retrieved objects.\n");

	std::shared_ptr<DwMaster> master(m->clone());
	DSPdebugMessage("Cloned the master.\n");

	/** create branching objects */
	std::shared_ptr<DspBranchObj> branch(new DspBranchObj);
	DSPdebugMessage("Created a branch object.\n");

	/** get reference point */
	std::vector<double> refsol;
	dynamic_cast<DwModel*>(model_)->getRefSol(refsol);

	// fix the solution
	fixSolution(tss, master.get(), branch.get(), refsol);

	master->setBranchingObjects(branch.get());
	DSPdebugMessage("Set a branch object.\n");

	int loglevel = master->getParPtr()->getIntParam("LOG_LEVEL");
	int dwevalub = master->getParPtr()->getIntParam("DW/EVAL_UB");
	double gaptol = master->getParPtr()->getDblParam("DW/GAPTOL");
	message->logLevel_ = 0;
	master->getParPtr()->setIntParam("DW/EVAL_UB", -1);
	master->getParPtr()->setDblParam("DW/GAPTOL", 0.0);

	DSP_RTN_CHECK_THROW(master->solve());
	DSPdebugMessage("Solved the master (status: %d).\n", master->getStatus());

	message->logLevel_ = loglevel;
	master->getParPtr()->setIntParam("DW/EVAL_UB", dwevalub);
	master->getParPtr()->setDblParam("DW/GAPTOL", gaptol);
	message->print(1, "Heuristic (Smip) returns %d.\n", master->getStatus());

	switch (master->getStatus()) {
	case DSP_STAT_OPTIMAL:
	case DSP_STAT_FEASIBLE:
	case DSP_STAT_LIM_ITERorTIME: {
		message->print(1, "Heuristic (Smip) found %e (vs. %e).\n", master->getBestDualObjective(), objective);
		if (master->getBestDualObjective() < objective) {
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

void DwSmip::fixSolution(TssModel* tss, DwMaster* master, DspBranchObj* branch, std::vector<double>& sol) {
	/** round and fix */
	for (int j = 0; j < tss->getNumCols(0); ++j) {
		if (master->ctype_orig_[j] != 'C') {
			/** round */
			double rounded = round(sol[j]);
			sol[j] = std::max(std::min(rounded, master->cubd_node_[j]), master->clbd_node_[j]);
		}
		/** fix */
		for (int s = 0; s < tss->getNumScenarios(); ++s)
			branch->push_back(tss->getNumCols(0) * s + j, sol[j], sol[j]);
	}
}
