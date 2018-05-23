/*
 * DwModel.cpp
 *
 *  Created on: Feb 13, 2017
 *      Author: kibaekkim
 */

//#define DSP_DEBUG

#include <DantzigWolfe/DwModel.h>
#include <DantzigWolfe/DwHeuristic.h>
#include <DantzigWolfe/DwBranchInt.h>
#include <DantzigWolfe/DwBranchNonant.h>
#include <Model/TssModel.h>

DwModel::DwModel(): DspModel(), branch_(NULL), infeasibility_(0.0) {}

DwModel::DwModel(DecSolver* solver): DspModel(solver), infeasibility_(0.0) {
	DwMaster* master = dynamic_cast<DwMaster*>(solver_);
	primsol_.resize(master->ncols_orig_);

	/** add heuristics */
	if (par_->getBoolParam("DW/HEURISTICS")) {
		if (par_->getBoolParam("DW/HEURISTICS/ROUNDING"))
			heuristics_.push_back(new DwRounding("Rounding", *this));
		if (par_->getBoolParam("DW/HEURISTICS/SMIP") && solver_->getModelPtr()->isStochastic())
			heuristics_.push_back(new DwSmip("Smip", *this));
	}

	switch (par_->getIntParam("DW/BRANCH")) {
		case BRANCH_INT:
			branch_ = new DwBranchInt(this);
			break;
		case BRANCH_NONANT:
			if (solver_->getModelPtr()->isStochastic())
				branch_ = new DwBranchNonant(this);
			else {
				solver_->getMessagePtr()->print(0, "Disabled nonanticipativity branching.\n");
				branch_ = new DwBranchInt(this);
			}
			break;
		default:
			branch_ = new DwBranchInt(this);
			break;
	}
}

DwModel::~DwModel() {
	for (unsigned i = 0; i < heuristics_.size(); ++i) {
		FREE_PTR(heuristics_[i]);
	}
	FREE_PTR(branch_);
}

DSP_RTN_CODE DwModel::solve() {
	BGN_TRY_CATCH

	DwMaster* master = dynamic_cast<DwMaster*>(solver_);
	DspMessage* message = master->getMessagePtr();

	/** set best primal objective value */
	master->setBestPrimalObjective(bestprimobj_);

	/** solve master */
	master->solve();

	status_ = master->getStatus();

	switch (status_) {
	case DSP_STAT_OPTIMAL:
	case DSP_STAT_FEASIBLE:
	case DSP_STAT_LIM_ITERorTIME: {

		primobj_ = master->getPrimalObjective();
		dualobj_ = master->getBestDualObjective();

		/** update best upper bound */
		if (master->getBestPrimalObjective() < bestprimobj_) {
			bestprimobj_ = master->getBestPrimalObjective();
			bestprimsol_.resize(master->ncols_orig_);
			for (int j = 0; j < master->ncols_orig_; ++j)
				bestprimsol_[j] = master->getBestPrimalSolutionOrig()[j];
			message->print(1, "Found new primal solution: %e\n", bestprimobj_);
		}

		if (primobj_ < 1.0e+20) {
			/** parse solution */
			int cpos = 0;
			std::fill(primsol_.begin(), primsol_.begin() + master->ncols_orig_, 0.0);
			for (auto it = master->cols_generated_.begin(); it != master->cols_generated_.end(); it++) {
				if ((*it)->active_) {
					if (fabs(master->getPrimalSolution()[cpos]) > 1.0e-10) {
						for (int i = 0; i < (*it)->x_.getNumElements(); ++i) {
							if ((*it)->x_.getIndices()[i] < master->ncols_orig_)
								primsol_[(*it)->x_.getIndices()[i]] += (*it)->x_.getElements()[i] * master->getPrimalSolution()[cpos];
						}
					}
					cpos++;
				}
			}
// #define WRITE_ALL_SOLS
#ifdef WRITE_ALL_SOLS
			{
				std::cout << "Number of columns: " << master->cols_generated_.size() << std::endl;
				std::ofstream fp_primsol("x.csv");
				for (auto it = master->cols_generated_.begin(); it != master->cols_generated_.end(); it++) {
					fp_primsol << (*it)->blockid_ << "," << (*it)->obj_;
					for (int j = 0; j < (*it)->x_.getNumElements(); ++j)
						fp_primsol << "," << (*it)->x_.getIndices()[j] << "," << (*it)->x_.getElements()[j];
					fp_primsol << std::endl;
				}
				fp_primsol.close();
			}
#endif
			//DspMessage::printArray(cpos, master->getPrimalSolution());

			/** calculate infeasibility */
			infeasibility_ = 0.0;
			for (int j = 0; j < master->ncols_orig_; ++j)
				if (master->ctype_orig_[j] != 'C') {
					infeasibility_ += fabs(primsol_[j] - floor(primsol_[j] + 0.5));
				}
			message->print(3, "Infeasibility: %+e\n", infeasibility_);

			bool isViolated = false;
			for (int j = 0; j < master->ncols_orig_; ++j) {
				double viol = std::max(master->clbd_node_[j] - primsol_[j], primsol_[j] - master->cubd_node_[j]);
				if (viol > 1.0e-6) {
					printf("Violated variable at %d by %e (%+e <= %+e <= %+e)\n", j, viol,
							master->clbd_node_[j], primsol_[j], master->cubd_node_[j]);
					isViolated = true;
				}
			}
			if (isViolated) throw "Invalid branching was performed.";

			/** run heuristics */
			if (par_->getBoolParam("DW/HEURISTICS") && infeasibility_ > 1.0e-6) {
				/** FIXME */
				for (auto it = heuristics_.begin(); it != heuristics_.end(); it++) {
					message->print(2, "Running [%s] heuristic:\n", (*it)->name());
					int found = (*it)->solution(bestprimobj_, bestprimsol_);
					if (found)
						message->print(2, "found bestprimobj %+e\n", bestprimobj_);
					else
						message->print(2, "Not found better primal solution\n");
				}
			}
		}

		break;
	}
	default:
		break;
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

bool DwModel::chooseBranchingObjects(std::vector<DspBranchObj*>& branchingObjs) {
	return branch_->chooseBranchingObjects(branchingObjs);
}
