//
// Created by Kibaek Kim on 8/27/16.
//

//#define DSP_DEBUG

/** standard */
#include <vector>
#include <algorithm>
/** coin */
#include "CoinTime.hpp"
/** Dsp */
#include "Utility/DspMessage.h"
#include "Model/BlkModel.h"

BlkModel::BlkModel() :
		nrows_full_(0),
		ncols_full_(0),
		nints_full_(0),
		primal_block_angular_(false),
		dual_block_angular_(false) {
	/** nothing to do */
}

BlkModel::BlkModel(const BlkModel& rhs) {
	Blocks blocks = const_cast<BlkModel&>(rhs).blocks();
	for (Blocks::iterator it = blocks.begin(); it != blocks.end(); ++it)
		addBlock(it->first, it->second);
}

DSP_RTN_CODE BlkModel::addBlock(int id, DetBlock* block) {
	/** check if id exists */
	if (blocks_.find(id) == blocks_.end())
		blocks_.insert(std::pair<const int,DetBlock*>(id, block));
	else {
		fprintf(stderr, "Block ID(%d) exists in the block model.\n", id);
		return DSP_RTN_ERR;
	}
	return DSP_RTN_OK;
}

DSP_RTN_CODE BlkModel::updateBlocks() {
#define FREE_MEMORY \
	FREE_ARRAY_PTR(ncols_coupled);

	int* ncols_coupled = NULL;

	BGN_TRY_CATCH

	std::vector<int> blockids;
	for (Blocks::iterator it = blocks_.begin(); it != blocks_.end(); ++it)
		blockids.push_back(it->first);

	double stime = CoinGetTimeOfDay();
	/** number of blocks */
	if (blockids.size() > 0) {
		/** coupling columns and rows */
		std::vector<int> master_coupling_cols;

		/** retrieve master information */
		DetBlock* master = block(blockids[0]);
		if (master == NULL)
			return DSP_RTN_ERR;
		const CoinPackedMatrix* master_mat = master->getConstraintMatrix();

		/** get number of columns, rows and integers */
		ncols_full_ = master->getNumCols();
		nrows_full_ = master->getNumRows();
		nints_full_ = master->getNumIntegers();

		/** to count how many subproblems are coupled with each column of the master. */
		ncols_coupled = new int [master->getNumCols()];
		CoinZeroN(ncols_coupled, master->getNumCols());

		for (unsigned i = 1; i < blockids.size(); ++i) {
			/** coupling columns and rows */
			std::vector<int> sub_coupling_cols;
			std::vector<int> sub_coupling_rows;

			/** retrieve subproblem */
			int id = blockids[i];
			DetBlock* sub = block(id);
			if (sub == NULL)
				return DSP_RTN_ERR;

			/** mark it is a dual block angular matrix */
			if (sub->getNumCols() > master->getNumCols())
				dual_block_angular_ = true;

			/** retrieve subproblem column indices */
			const CoinPackedMatrix* sub_mat = sub->getConstraintMatrix();
			int sub_num_col_indices = sub_mat->isColOrdered() ?
					sub_mat->getNumCols() : sub_mat->getNumElements();
			int* sub_col_indices = new int [sub_num_col_indices];
			if (sub_mat->isColOrdered())
				CoinCopyN(sub_mat->getVectorStarts(), sub_num_col_indices, sub_col_indices);
			else
				CoinCopyN(sub_mat->getIndices(), sub_num_col_indices, sub_col_indices);

			/** add number of columns, rows and integers */
			ncols_full_ += sub->getNumCols() - master->getNumCols();
			nrows_full_ += sub->getNumRows() - master->getNumRows();
			nints_full_ += sub->getNumIntegers() - master->getNumIntegers();

			/** sort indices to speed up the search */
			std::sort(sub_col_indices, sub_col_indices + sub_num_col_indices);

			/** find coupling columns and rows */
			for (int j = 0; j < sub_num_col_indices; ++j) {
				if (j > 0 && sub_col_indices[j-1] == sub_col_indices[j])
					continue;
				if (sub_col_indices[j] < master_mat->getNumCols()) {
					master_coupling_cols.push_back(sub_col_indices[j]);
					sub_coupling_cols.push_back(sub_col_indices[j]);
					/** find all the coupling rows */
					for (int row = 0; row < master->getNumRows(); ++row) {
						if (fabs(master_mat->getCoefficient(row, sub_col_indices[j])) > 1.0e-8)
							sub_coupling_rows.push_back(row);
					}
				}
			}

			/** erase duplicates */
			std::sort(sub_coupling_cols.begin(), sub_coupling_cols.end());
			sub_coupling_cols.erase(
					std::unique(sub_coupling_cols.begin(), sub_coupling_cols.end()),
					sub_coupling_cols.end());
			std::sort(sub_coupling_rows.begin(), sub_coupling_rows.end());
			sub_coupling_rows.erase(
					std::unique(sub_coupling_rows.begin(), sub_coupling_rows.end()),
					sub_coupling_rows.end());

			/** set coupling columns and rows */
			sub->setCouplingCols(sub_coupling_cols.size(), &sub_coupling_cols[0]);
			sub->setCouplingRows(sub_coupling_rows.size(), &sub_coupling_rows[0]);
			DSPdebugMessage("Coupling columns of block %d:\n", id);
			DSPdebug(DspMessage::printArray(sub->getNumCouplingCols(), sub->getCouplingCols()));
			DSPdebugMessage("Coupling rows of block %d:\n", id);
			DSPdebug(DspMessage::printArray(sub->getNumCouplingRows(), sub->getCouplingRows()));

			/** count the master columns coupled */
			for (unsigned j = 0; j < sub_coupling_cols.size(); ++j)
				ncols_coupled[sub_coupling_cols[j]]++;
		}

		/** erase duplicates */
		std::sort(master_coupling_cols.begin(), master_coupling_cols.end());
		master_coupling_cols.erase(
				std::unique(master_coupling_cols.begin(), master_coupling_cols.end()),
				master_coupling_cols.end());

		/** set coupling columns */
		master->setCouplingCols(master_coupling_cols.size(), &master_coupling_cols[0]);
		/** set coupling rows */
		int* master_coupling_rows = new int [master->getNumRows()];
		for (int i = 0; i < master->getNumRows(); ++i)
			master_coupling_rows[i] = i;
		master->setCouplingRows(master->getNumRows(), master_coupling_rows);
		FREE_ARRAY_PTR(master_coupling_rows);

		/** check if the full matrix is of a primal block angular form. */
		primal_block_angular_ = true;
		for (int j = 0; j < master->getNumCols(); ++j)
			if (ncols_coupled[j] > 1)
				primal_block_angular_ = false;

		if (primal_block_angular_ == true)
			printf("The constraint matrix is of a primal block angular.\n");
		if (dual_block_angular_ == true)
			printf("The constraint matrix is of a dual block angular.\n");
	}
	DSPdebugMessage("Update block time: %.4f\n", CoinGetTimeOfDay() - stime);

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

DetBlock* BlkModel::block(int id) {
	DetBlock* b = NULL;
	/** check if id exists */
	if (blocks_.find(id) == blocks_.end())
		fprintf(stderr, "Block ID(%d) does not exist in the block model.\n", id);
	else
		b = blocks_[id];
	return b;
}
