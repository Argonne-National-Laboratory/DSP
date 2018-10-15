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

	for (int j = 0; j < tss_->getNumCols(0); ++j) {
		devsol[j] = CoinMin(fabs(refsol[j] - master_->clbd_node_[j]), fabs(refsol[j] - master_->cubd_node_[j]));
		DSPdebugMessage("col %d: devsol %e\n", j, devsol[j]);
	}

#if 0
	std::vector<double> diffsol(tss_->getNumCols(0), 0.0);
	/** use l2-norm */
	for (auto it = master_->cols_generated_.begin(); it != master_->cols_generated_.end(); it++) {
		if ((*it)->active_) {
			int master_index = (*it)->master_index_;
			int sind = (*it)->blockid_;
			double weight = master_->getPrimalSolution()[master_index];
			if (weight > 1.0e-10) {
				for (int i = 0; i < (*it)->x_.getNumElements(); ++i) {
					if ((*it)->x_.getIndices()[i] < tss_->getNumCols(0) * tss_->getNumScenarios()) {
						int k = (*it)->x_.getIndices()[i] % tss_->getNumCols(0);
						diffsol[k] += pow((*it)->x_.getElements()[i] - refsol[k], 2.0) * weight * tss_->getProbability()[sind];
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
#endif
}
