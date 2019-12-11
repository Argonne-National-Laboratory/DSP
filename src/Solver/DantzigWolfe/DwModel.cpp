/*
 * DwModel.cpp
 *
 *  Created on: Feb 13, 2017
 *      Author: kibaekkim
 */

//#define DSP_DEBUG

#include "Solver/DantzigWolfe/DwModel.h"
#include "Solver/DantzigWolfe/DwHeuristic.h"
#include "Solver/DantzigWolfe/DwBranchInt.h"
#include "Solver/DantzigWolfe/DwBranchNonant.h"
#include "Solver/DantzigWolfe/DwBranchNonant2.h"
#include "Solver/DantzigWolfe/DwBranchGenDisj.h"
#include "Model/TssModel.h"

DwModel::DwModel(): DspModel(), branch_(NULL) {}

DwModel::DwModel(DecSolver* solver): DspModel(solver) {
	DwMaster* master = dynamic_cast<DwMaster*>(solver_);
	primsol_.resize(master->ncols_orig_);

	/** add heuristics */
	if (par_->getBoolParam("DW/HEURISTICS")) {
		if (par_->getBoolParam("DW/HEURISTICS/ROUNDING"))
			heuristics_.push_back(new DwRounding("Rounding", *this));
		if (solver_->getModelPtr()->isStochastic() == false)
			par_->setBoolParam("DW/HEURISTICS/SMIP", false);
		if (par_->getBoolParam("DW/HEURISTICS/SMIP"))
			heuristics_.push_back(new DwSmip("Smip", *this));
	}

	if (par_->getBoolParam("DW/SPLIT_VARS")) {
		par_->setIntParam("DW/BRANCH", BRANCH_NONANT);
		branch_ = new DwBranchNonant(this);
	} else {
		if (solver_->getModelPtr()->isStochastic()) {
			switch (par_->getIntParam("DW/BRANCH")) {
				case BRANCH_NONANT:
					branch_ = new DwBranchNonant(this);
					break;
				case BRANCH_NONANT2:
					branch_ = new DwBranchNonant2(this);
					break;
				case BRANCH_DISJUNCTION_TEST:
					printf("\n*** WARNING: This is a testing code for branching on general disjunctions. ***\n\n");
					branch_ = new DwBranchGenDisj(this);
					break;
				case BRANCH_INT:
				default:
					branch_ = new DwBranchInt(this);
					break;
			}
		} else {
			if (par_->getIntParam("DW/BRANCH") != BRANCH_INT)
				printf("\n*** WARNING: Changed the branching method. The problem structure supports the integer branching method only. ***\n\n");
			branch_ = new DwBranchInt(this);
		}
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

		/** parse solution */
		DSP_RTN_CHECK_RTN_CODE(parsePrimSolution());

		/** run heuristics */
		if (par_->getBoolParam("DW/HEURISTICS") && infeasibility_ > 1.0e-6) {
			/** FIXME */
			std::vector<double> primsol(primsol_);
			for (auto it = heuristics_.begin(); it != heuristics_.end(); it++) {
				primsol_ = primsol;
				message->print(2, "Running [%s] heuristic:\n", (*it)->name());
				int found = (*it)->solution(bestprimobj_, bestprimsol_);
				if (found)
					message->print(2, "found bestprimobj %+e\n", bestprimobj_);
				else
					message->print(2, "Not found better primal solution\n");
			}
			primsol_ = primsol;
		}

		break;
	}
	case DSP_STAT_LIM_DUAL_OBJ:
		dualobj_ = master->getBestDualObjective();
		break;
	default:
		break;
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DwModel::parsePrimSolution() {
	BGN_TRY_CATCH

	DwMaster* master = dynamic_cast<DwMaster*>(solver_);
	DspMessage* message = master->getMessagePtr();

	if (par_->getBoolParam("DW/SPLIT_VARS")) {
		/** NOTE: 
			This is purely for reporting numerical experiment to paper.
			Using the restricted master solution makes more sense in practice. */
		DSP_RTN_CHECK_RTN_CODE(parseLastIterSolution());

		/** set infinite infeasibility */
		infeasibility_ = COIN_DBL_MAX;
	} else {
		if (solver_->getModelPtr()->isStochastic() == true && 
			par_->getIntParam("DW/BRANCH") == BRANCH_NONANT) {
			/** NOTE: 
				This is purely for reporting numerical experiment to paper.
				Using the restricted master solution makes more sense in practice. */
			DSP_RTN_CHECK_RTN_CODE(parseLastIterSolution());

			/** set infinite infeasibility */
			infeasibility_ = COIN_DBL_MAX;
		} else {
			if (primobj_ < 1.0e+20) {
				/** getting a Dantzig-Wolfe solution */
				DSP_RTN_CHECK_RTN_CODE(parseDantzigWolfeSolution());

				/** calculate infeasibility */
				infeasibility_ = 0.0;
				for (int j = 0; j < master->ncols_orig_; ++j)
					if (master->ctype_orig_[j] != 'C')
						infeasibility_ += fabs(primsol_[j] - floor(primsol_[j] + 0.5));
				message->print(1, "Infeasibility: %+e\n", infeasibility_);

				/** extra checking for column bounds */
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
			} else {
				/** disable running heuristics */
				infeasibility_ = 0.0;
			}
		}

	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

DSP_RTN_CODE DwModel::parseDantzigWolfeSolution() {
	BGN_TRY_CATCH

	DwMaster* master = dynamic_cast<DwMaster*>(solver_);

	std::fill(primsol_.begin(), primsol_.begin() + master->ncols_orig_, 0.0);

	for (auto it = master->cols_generated_.begin(); it != master->cols_generated_.end(); it++) {
		if ((*it)->active_) {
			int master_index = (*it)->master_index_;
			if (fabs(master->getPrimalSolution()[master_index]) > 1.0e-10) {
				for (int i = 0; i < (*it)->x_.getNumElements(); ++i) {
					if ((*it)->x_.getIndices()[i] < master->ncols_orig_)
						primsol_[(*it)->x_.getIndices()[i]] += (*it)->x_.getElements()[i] * master->getPrimalSolution()[master_index];
				}
			}
		}
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

DSP_RTN_CODE DwModel::parseLastIterSolution() {
	BGN_TRY_CATCH

	DwMaster* master = dynamic_cast<DwMaster*>(solver_);

	std::fill(primsol_.begin(), primsol_.begin() + master->ncols_orig_, 0.0);

	for (unsigned s = 0; s < master->getLastSubprobSolutions().size(); ++s) {
		const CoinPackedVector* sol = master->getLastSubprobSolutions()[s];
		int sNumElements = sol->getNumElements();
		const int* sIndices = sol->getIndices();
		const double* sElements = sol->getElements();
		for (int j = 0; j < sNumElements; ++j)
			primsol_[sIndices[j]] += sElements[j];
	}

	/** TODO: we can calculuate infeasibility here with respect to the nonanticipativity */

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

bool DwModel::chooseBranchingObjects(std::vector<DspBranchObj*>& branchingObjs) {
	return branch_->chooseBranchingObjects(branchingObjs);
}
