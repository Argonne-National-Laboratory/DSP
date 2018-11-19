/*
 * DwBranchNonant.cpp
 *
 *  Created on: May 2, 2018
 *      Author: Kibaek Kim
 */

//#define DSP_DEBUG
#include <DantzigWolfe/DwBranchNonant.h>

DwBranchNonant::DwBranchNonant(DwModel* model) : DwBranch(model) {
	DecSolver* solver = model_->getSolver();
	master_ = dynamic_cast<DwMaster*>(solver);
	if (solver->getModelPtr()->isStochastic())
		tss_ = dynamic_cast<TssModel*>(solver->getModelPtr());
}

bool DwBranchNonant::chooseBranchingObjects(
		std::vector<DspBranchObj*>& branchingObjs /**< [out] branching objects */) {

	if (tss_ == NULL) return false;

	bool branched = false;
	int branchingIndex = -1;
	double branchingValue, branchingDownValue, branchingUpValue;

	std::vector<double> refsol; // reference solution
	std::vector<double> devsol; // devations from the refsol
	std::vector<std::vector<double>> densesol; // dense solution vector

	DspBranchObj* branchingUp = NULL;
	DspBranchObj* branchingDn = NULL;

	BGN_TRY_CATCH

	/** retreive message object */
	DspMessage* message = model_->getSolver()->getMessagePtr();

	/** calculate reference solution */
	getRefSol(refsol);

	/** calculate variances from the reference solution */
	getDevSol(refsol, devsol);

#ifdef DSP_DEBUG
	printf("Reference Solution:\n");
	DspMessage::printArray(refsol.size(), &refsol[0]);
	printf("Deviations:\n");
	DspMessage::printArray(devsol.size(), &devsol[0]);
#endif

	double maxdev = 0.0;
	for (int j = 0; j < tss_->getNumCols(0); ++j) {
		if (devsol[j] > CoinMax(epsilon_, maxdev)) {
			maxdev = devsol[j];
			branchingIndex = j;
			branchingValue = refsol[j];
			if (tss_->getCtypeCore(0)[j] == 'C') {
				branchingDownValue = refsol[j] - epsilonBB_;
				branchingUpValue = refsol[j] + epsilonBB_;
			} else {
				branchingDownValue = floor(refsol[j]);
				branchingUpValue = ceil(refsol[j]);
			}
		}
	}
	DSPdebugMessage("maxdev %e\n", maxdev);

	if (branchingIndex > -1) {

		branched = true;

		message->print(2, "Creating branch objects on column %d (value %e): [%e,%e] and [%e,%e]\n",
			branchingIndex, branchingValue, master_->clbd_node_[branchingIndex], branchingDownValue, branchingUpValue, master_->cubd_node_[branchingIndex]);

		/** creating branching objects */
		branchingUp = new DspBranchObj();
		branchingDn = new DspBranchObj();
		for (int j = 0; j < tss_->getNumCols(0) * tss_->getNumScenarios(); ++j) {
			if (branchingIndex == j % tss_->getNumCols(0)) {
				branchingUp->push_back(j, CoinMin(branchingUpValue, master_->cubd_node_[j]), master_->cubd_node_[j]);
				branchingDn->push_back(j, master_->clbd_node_[j], CoinMax(master_->clbd_node_[j], branchingDownValue));
			} else if (master_->clbd_node_[j] > master_->clbd_orig_[j] || master_->cubd_node_[j] < master_->cubd_orig_[j]) {
				/** store any bound changes made in parent nodes */
				//DSPdebugMessage("Adjusting bound change on column %d: [%e,%e]\n", j, master_->clbd_node_[j], master_->cubd_node_[j]);
				branchingUp->push_back(j, master_->clbd_node_[j], master_->cubd_node_[j]);
				branchingDn->push_back(j, master_->clbd_node_[j], master_->cubd_node_[j]);
			}
		}

		/** set best dual bounds */
		branchingUp->bestBound_ = master_->getBestDualObjective();
		branchingDn->bestBound_ = master_->getBestDualObjective();

		/** assign best dual solutions */
		branchingUp->dualsol_.assign(master_->getBestDualSolution(), master_->getBestDualSolution() + master_->nrows_);
		branchingDn->dualsol_.assign(master_->getBestDualSolution(), master_->getBestDualSolution() + master_->nrows_);

		/** set branching directions */
		branchingUp->direction_ = 1;
		branchingDn->direction_ = -1;

		branchingUp->solEstimate_ = maxdev;
		branchingDn->solEstimate_ = maxdev;

		/** add branching objects */
		branchingObjs.push_back(branchingUp);
		branchingObjs.push_back(branchingDn);
	} else {
		DSPdebugMessage("No branch object is found.\n");
	}

	END_TRY_CATCH_RTN(;,false)

	return branched;
}

void DwBranchNonant::getRefSol(std::vector<double>& refsol) {
	if (model_->getParPtr()->getBoolParam("DW/SPLIT_VARS") &&
		model_->getSolver()->getModelPtr()->isStochastic() == false) {
		/** TODO: Liu's implementation goes here. */
		/** TODO: This is the place for implementing heuristics for energey storage application
			Algorithmic Steps
				1. access the master block
				2. collect coupling columns
				3. collect the solutions for coupling columns
				4. weight average
		*/
		DecBlkModel* dec = dynamic_cast<DecBlkModel*>(model_);
		DetBlock* master_block = dec->blkPtr()->block(0);
		const CoinPackedMatrix* master_mat = master_block->getConstraintMatrix();

		/** set weight parameter for time-related variables in coupling constraints */
		double weight_con_cols = 0.5;

		/** collect coupling columns with nonzero coefficient */
		std::vector<int> coupling_cols(master_mat->getNumElements(), 0);
		/** collect coupling columns with nonzero coefficient with duplicates*/
		std::vector<int> coupling_cols_dup(master_mat->getNumElements(), 0);
		/** get all column indices from the master matrix */
		for (int i = 0; i < master_mat->getNumElements(); ++i){
			coupling_cols[i] = master_mat->getIndices()[i];
			coupling_cols_dup[i] = master_mat->getIndices()[i];
		}
		/** remove duplicates */
		std::sort(coupling_cols.begin(), coupling_cols.end());
		coupling_cols.erase(std::unique(coupling_cols.begin(), coupling_cols.end()), coupling_cols.end());
		printf("master block: total columns %d coupling columns: %d\n", master_block->getNumCols(), coupling_cols.size());
		/** collect coupling solutions */
		std::vector<double> coupling_cols_solution(coupling_cols.size(),0.0);
		for (int i=0; i < coupling_cols.size(); ++i) {
			coupling_cols_solution[i] = model_->getPrimalSolution()[coupling_cols[i]];
		}
		/** weight average */
		/** collect the equalities */
		std::vector<std::vector <int>> coupling_equal
		coupling_equal[0].push_back(coupling_cols_dup[0]);
		coupling_equal[0].push_back(coupling_cols_dup[1]);
		int count_num_equal = 0;
		coupling_cols_dup.erase(coupling_cols_dup.begin(),coupling_cols_dup.begin()+1);
		while (coupling_cols_dup.size()>0) {
			/**double check if to refer the last element of vector is vec.end*/
			const int* pos = std::find(coupling_cols_dup.begin(), coupling_cols_dup.end(), coupling_equal[count_num_equal].back());
			if (pos != coupling_cols_dup.end()) {
				coupling_equal[count_num_equal].push_back(coupling_cols_dup[pos]);
				coupling_cols_dup.erase(coupling_cols_dup.begin()+pos-1,coupling_cols_dup.begin()+pos);
			}
			else {
				++count_num_equal;
				coupling_equal[count_num_equal].push_back(coupling_cols_dup[0]);
				coupling_equal[count_num_equal].push_back(coupling_cols_dup[1]);
				coupling_cols_dup.erase(coupling_cols_dup.begin(),coupling_cols_dup.begin()+1);
			}

		}
		/** calculate and assign the average values of equal variables*/
		refsol.resize(coupling_equal.size(),0.0);
		for (int i = 0; i < coupling_equal.size(); ++i) {
			double sum_temp = 0;
			for (int j = 0; j < coupling_equal[i].size(); ++j) {
				auto pos_coup = std::find(coupling_cols.begin(), coupling_cols.end(), coupling_equal[i][j]);
				sum_temp += coupling_cols_solution[distance(coupling_cols.begin(), pos_coup);];
			}
			sum_temp = sum_temp/double(coupling_equal.size());

			/** assign averaged value to variables */
			for (int j = 0; j < coupling_equal[i].size(); ++j) {
				auto pos_coup = std::find(coupling_cols.begin(), coupling_cols.end(), coupling_equal[i][j]);
				refsol[distance(coupling_cols.begin(),pos_coup)] = sum_temp;
			}
		}

	} else {
		refsol.resize(tss_->getNumCols(0), 0.0);
		for (int j = 0; j < tss_->getNumCols(0) * tss_->getNumScenarios(); ++j) {
			int s = j / tss_->getNumCols(0);
			refsol[j % tss_->getNumCols(0)] += model_->getPrimalSolution()[j] * tss_->getProbability()[s];
		}
	}
}

void DwBranchNonant::getDevSol(std::vector<double>& refsol, std::vector<double>& devsol) {
	if (model_->getParPtr()->getBoolParam("DW/SPLIT_VARS") &&
		model_->getSolver()->getModelPtr()->isStochastic() == false) {
		/** TODO: Liu's implementation goes here. */
		devsol.resize(coupling_equal.size(),0.0);
		std::vector<double> diffsol(coupling_equal.size(), 0.0);
	} else {
		devsol.resize(tss_->getNumCols(0), 0.0);
//#define USE_TWONORM
#ifdef USE_TWONORM
		std::vector<double> diffsol(tss_->getNumCols(0), 0.0);
		/** use l2-norm */
		for (unsigned s = 0; s < master_->getLastSubprobSolutions().size(); ++s) {
			const CoinPackedVector* sol = master_->getLastSubprobSolutions()[s];
			int sNumElements = sol->getNumElements();
			const int* sIndices = sol->getIndices();
			const double* sElements = sol->getElements();
			for (int j = 0; j < sNumElements; ++j)
				if (sIndices[j] < tss_->getNumCols(0) * tss_->getNumScenarios()) {
					int k = sIndices[j] % tss_->getNumCols(0);
					diffsol[k] += pow(sElements[j] - refsol[k], 2.0) * tss_->getProbability()[s];
				}
		}
		for (int k = 0; k < tss_->getNumCols(0); ++k)
			devsol[k] = diffsol[k] > 1.0e-10 ? sqrt(diffsol[k]) : 0.0;
#else
		std::vector<double> maxsol(tss_->getNumCols(0), -COIN_DBL_MAX);
		std::vector<double> minsol(tss_->getNumCols(0), +COIN_DBL_MAX);
		/** calculate max value first */
		for (unsigned s = 0; s < master_->getLastSubprobSolutions().size(); ++s) {
			const CoinPackedVector* sol = master_->getLastSubprobSolutions()[s];
			int sNumElements = sol->getNumElements();
			const int* sIndices = sol->getIndices();
			const double* sElements = sol->getElements();
			for (int j = 0; j < sNumElements; ++j)
				if (sIndices[j] < tss_->getNumCols(0) * tss_->getNumScenarios()) {
					int k = sIndices[j] % tss_->getNumCols(0);
					maxsol[k] = CoinMax(maxsol[k], sElements[j]);
					minsol[k] = CoinMin(minsol[k], sElements[j]);
				}
		}
		for (int k = 0; k < tss_->getNumCols(0); ++k)
			devsol[k] = CoinMax(maxsol[k] - minsol[k], 0.0);
#endif
	}
}

#undef USE_TWONORM
