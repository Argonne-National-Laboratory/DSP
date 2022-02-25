/*
 * DwModelSmip.cpp
 *
 *  Created on: Jan 20, 2020
 *      Author: Kibaek Kim
 */

//#define DSP_DEBUG

#include "Solver/DantzigWolfe/DwModelSmip.h"
#include "Solver/DantzigWolfe/DwHeuristic.h"
#include "Solver/DantzigWolfe/DwBranchInt.h"
#include "Solver/DantzigWolfe/DwBranchNonant.h"
#include "Solver/DantzigWolfe/DwBranchNonant2.h"
#include "Solver/DantzigWolfe/DwBranchGenDisj.h"

DwModelSmip::DwModelSmip():
DwModel() {
    tss_ = dynamic_cast<TssModel*>(solver_->getModelPtr());
}

DwModelSmip::DwModelSmip(DecSolver* solver):
DwModel(solver) {
    tss_ = dynamic_cast<TssModel*>(solver_->getModelPtr());
	// printf("DwModelSmip constructor.\n");
}

DwModelSmip::~DwModelSmip() {
    tss_ = NULL;
}

DSP_RTN_CODE DwModelSmip::initBranch() {
	BGN_TRY_CATCH

    switch (par_->getIntParam("DW/BRANCH")) {
        case BRANCH_NONANT:
            branch_ = new DwBranchNonant(this);
            break;
        case BRANCH_NONANT2:
			// printf("Created BRANCH_NONANT2 rule.\n");
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

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DwModelSmip::initHeuristic() {
	BGN_TRY_CATCH

	if (par_->getBoolParam("DW/HEURISTICS")) {
		if (par_->getBoolParam("DW/HEURISTICS/ROUNDING"))
			heuristics_.push_back(new DwRounding("Rounding", *this));
		if (par_->getBoolParam("DW/HEURISTICS/SMIP"))
			heuristics_.push_back(new DwSmip("Smip", *this));
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DwModelSmip::parsePrimSolution() {
	BGN_TRY_CATCH

	DwMaster* master = dynamic_cast<DwMaster*>(solver_);
	DspMessage* message = master->getMessagePtr();

	if (par_->getIntParam("DW/BRANCH") == BRANCH_NONANT) {
		/** NOTE: 
			This is purely for reporting numerical experiment to paper.
			Using the restricted master solution makes more sense in practice. */
		DSP_RTN_CHECK_RTN_CODE(parseLastIterSolution());
	} else {
        DwModel::parsePrimSolution();

        /** NOTE: This may not be necessary. */
        if (par_->getIntParam("DW/BRANCH") == BRANCH_NONANT2 &&
            infeasibility_ > feastol_) {

            std::vector<double> primsol(master->ncols_orig_, 0.0);

            // find the index for the maximum alpha
            std::vector<int> maxindex(tss_->getNumScenarios(), -1);
            std::vector<double> maxalpha(tss_->getNumScenarios(), 0.0);
            for (unsigned i = 0; i < master->cols_generated_.size(); ++i) {
                DwCol* col = master->cols_generated_[i];
                if (col->active_) {
                    int blockid = col->blockid_;
                    int master_index = col->master_index_;
                    if (master->getPrimalSolution()[master_index] > maxalpha[blockid]) {
                        maxindex[blockid] = i;
                        maxalpha[blockid] = master->getPrimalSolution()[master_index];
                    }
                }
            }

            for (int s = 0; s < tss_->getNumScenarios(); ++s) {
                DwCol* col = master->cols_generated_[maxindex[s]];
                for (int i = 0; i < col->x_.getNumElements(); ++i) {
                    if (col->x_.getIndices()[i] < master->ncols_orig_)
                        primsol[col->x_.getIndices()[i]] = col->x_.getElements()[i];
                }
            }

			// printf("First-stage solution:\n");
			// DspMessage::printArray(tss_->getNumCols(0) * tss_->getNumScenarios(), &primsol[0]);

            // calculate infeasibility_
            double infeasibility = getMaxDev(primsol);
            message->print(1, "Infeasibility (Nonanti): %+e\n", infeasibility);

            if (infeasibility <= feastol_) {
                primsol_ = primsol;
				infeasibility_ = infeasibility;
            }
        }
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

DSP_RTN_CODE DwModelSmip::parseLastIterSolution() {
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
	// calculate infeasibility_
	infeasibility_ = getMaxDev(primsol_);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

void DwModelSmip::getRefSol(std::vector<double>& refsol) {
	if (par_->getIntParam("DW/BRANCH") == BRANCH_NONANT) {
		std::fill(refsol.begin(), refsol.begin() + tss_->getNumCols(0) * tss_->getNumScenarios(), 0.0);
		for (int j = 0; j < tss_->getNumCols(0) * tss_->getNumScenarios(); ++j) {
			int s = j / tss_->getNumCols(0);
			refsol[j % tss_->getNumCols(0)] += primsol_[j] * tss_->getProbability()[s];
		}
		for (int s = 1; s < tss_->getNumScenarios(); ++s) {
			std::copy(refsol.begin(), refsol.begin() + tss_->getNumCols(0), 
				refsol.begin() + tss_->getNumCols(0) * s);
		}
	} else {
        refsol = primsol_;
    }
}

double DwModelSmip::getMaxDev(std::vector<double> primsol) {
	double maxdev = 0.0;

	std::vector<double> maxsol(tss_->getNumCols(0), -COIN_DBL_MAX);
	std::vector<double> minsol(tss_->getNumCols(0), +COIN_DBL_MAX);

	/** calculate max value first */
	for (int j = 0; j < tss_->getNumCols(0) * tss_->getNumScenarios(); ++j) {
		int k = j % tss_->getNumCols(0);
		maxsol[k] = CoinMax(maxsol[k], primsol[j]);
		minsol[k] = CoinMin(minsol[k], primsol[j]);
	}
	for (int k = 0; k < tss_->getNumCols(0); ++k) {
		DSPdebugMessage("k %d maxsol %e minsol %e\n", k, maxsol[k], minsol[k]);
		maxdev = CoinMax(maxsol[k] - minsol[k], maxdev);
	}

	return maxdev;
}