/*
 * DwHeuristic.cpp
 *
 *  Created on: Feb 12, 2017
 *      Author: kibaekkim
 */

#include <DantzigWolfe/DwHeuristic.h>
#include "Solver/DantzigWolfe/DwMaster.h"

DwRounding::DwRounding(std::string name, DecSolver &solver):
DspHeuristic(name, solver) {}

int DwRounding::solution(double &objective, std::vector<double> &solution) {

	int found = 0;
	std::vector<int> rowadded;
	std::vector<int> rowind;
	std::vector<double> rowval;
	double rounded;

	DwMaster* master = dynamic_cast<DwMaster*>(solver_);
	OsiSolverInterface* si = master->getSiPtr();

	double primobj = master->getPrimalObjective();
	double dualobj = master->getDualObjective();
	int status = master->getStatus();
	std::vector<double> primsol(master->getPrimalSolution(),
			master->getPrimalSolution() + master->ncols_orig_);

	master->setPrimalObjective(COIN_DBL_MAX);

	for (int j = 0; j < master->ncols_orig_; ++j) {
		if (master->org_ctype_[j] != 'C') {
			/** round */
			rounded = round(primsol[j]);
			rounded = std::min(rounded, master->node_cubd_[j]);
			rounded = std::max(rounded, master->node_clbd_[j]);
			/** fix */
			rowind.clear();
			rowval.clear();
			int pos = 0;
			for (auto it = master->cols_generated_.begin();
					it != master->cols_generated_.end(); it++) {
				if ((*it)->active_) {
					int sparse_index = (*it)->x_.findIndex(j);
					if (sparse_index > -1) {
						rowind.push_back(pos);
						rowval.push_back((*it)->x_.getElements()[sparse_index]);
					}
					pos++;
				}
			}
			rowadded.push_back(si->getNumRows());
			si->addRow(rowind.size(), &rowind[0], &rowval[0], rounded, rounded);
		}
	}

	DSP_RTN_CHECK_THROW(master->solve());

	switch (master->getStatus()) {
	case DSP_STAT_OPTIMAL:
	case DSP_STAT_FEASIBLE:
	case DSP_STAT_LIM_ITERorTIME:
		if (master->getPrimalObjective() < objective) {
			objective = master->getPrimalObjective();
			solution.assign(master->getPrimalSolution(), master->getPrimalSolution() + master->ncols_orig_);
			found = 1;
		}
		break;
	default:
		break;
	}

	/** switch to phase 2 */
	DSP_RTN_CHECK_RTN_CODE(master->switchToPhase2());

	/** delete rows added */
	si->deleteRows(rowadded.size(), &rowadded[0]);

	/** restore original problem data */
	master->setBestPrimalObjective(objective);
	master->setPrimalObjective(primobj);
	master->setDualObjective(dualobj);
	master->setPrimalSolution(&primsol[0]);
	master->setStatus(status);

	return found;
}
