/*
 * DwHeuristic.cpp
 *
 *  Created on: Feb 12, 2017
 *      Author: kibaekkim
 */

#include <memory>
#include <DantzigWolfe/DwHeuristic.h>
#include <DantzigWolfe/DwMaster.h>

DwRounding::DwRounding(std::string name, DecSolver &solver):
DspHeuristic(name, solver) {}

int DwRounding::solution(double &objective, std::vector<double> &solution) {

	/** create branching objects */
	std::shared_ptr<DspBranch> branch(new DspBranch);

	int found = 0;
	double rounded;

	DwMaster* master = dynamic_cast<DwMaster*>(solver_);
	OsiSolverInterface* si = master->getSiPtr();

	double primobj = master->getPrimalObjective();
	double dualobj = master->getDualObjective();
	int status = master->getStatus();
	std::vector<double> primsol(master->getPrimalSolution(), master->getPrimalSolution() + master->ncols_orig_);

	master->setPrimalObjective(COIN_DBL_MAX);

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

	DSP_RTN_CHECK_THROW(master->solve());

	switch (master->getStatus()) {
	case DSP_STAT_OPTIMAL:
	case DSP_STAT_FEASIBLE:
	case DSP_STAT_LIM_ITERorTIME:
		if (master->getPrimalObjective() < objective) {
			objective = master->getBestDualObjective();
			solution.assign(master->getPrimalSolution(), master->getPrimalSolution() + master->ncols_orig_);
			found = 1;
		}
		break;
	default:
		break;
	}

	/** switch to phase 2 */
	DSP_RTN_CHECK_RTN_CODE(master->switchToPhase2());

	/** restore original problem data */
	master->setBestPrimalObjective(objective);
	master->setPrimalObjective(primobj);
	master->setDualObjective(dualobj);
	master->setPrimalSolution(&primsol[0]);
	master->setStatus(status);

	return found;
}
