/*
 * DwBranchNonant.cpp
 *
 *  Created on: May 8, 2018
 *      Author: Kibaek Kim
 */

//#define DSP_DEBUG
#include <DantzigWolfe/DwBranchNonant2.h>
#include <Model/TssModel.h>

void DwBranchNonant2::getRefSol(std::vector<double>& refsol) {
	refsol = model_->getPrimalSolution();
}

void DwBranchNonant2::getDevSol(std::vector<double>& refsol, std::vector<double>& devsol) {
	devsol.resize(tss_->getNumCols(0), 0.0);
	//std::vector<double> maxsol(tss_->getNumCols(0), -COIN_DBL_MAX);
	//std::vector<double> minsol(tss_->getNumCols(0), +COIN_DBL_MAX);
	std::vector<double> diffsol(tss_->getNumCols(0), 0.0);

	/** use l2-norm */
	for (auto it = master_->cols_generated_.begin(); it != master_->cols_generated_.end(); it++) {
		if ((*it)->active_) {
			int master_index = (*it)->master_index_;
			double weight = master_->getPrimalSolution()[master_index];
			if (weight > 1.0e-10) {
				for (int i = 0; i < (*it)->x_.getNumElements(); ++i) {
					if ((*it)->x_.getIndices()[i] < tss_->getNumCols(0) * tss_->getNumScenarios()) {
						int k = (*it)->x_.getIndices()[i] % tss_->getNumCols(0);
						//maxsol[k] = CoinMax(maxsol[k], (*it)->x_.getElements()[i]);
						//minsol[k] = CoinMin(minsol[k], (*it)->x_.getElements()[i]);
						diffsol[k] += pow((*it)->x_.getElements()[i] - refsol[k], 2.0) * weight;
					}
				}
			}
		}
	}
	for (int k = 0; k < tss_->getNumCols(0); ++k) {
		//devsol[k] = CoinMax(maxsol[k] - minsol[k], 0.0);
		//DSPdebugMessage("col %d: devsol %e maxsol %e minsol %e\n", k, devsol[k], maxsol[k], minsol[k]);
		devsol[k] = diffsol[k] > 1.0e-10 ? sqrt(diffsol[k]) : 0.0;
		DSPdebugMessage("col %d: devsol %e diffsol %e\n", k, devsol[k], diffsol[k]);
	}
}
