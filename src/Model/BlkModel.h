//
// Created by Kibaek Kim on 8/27/16.
//

#ifndef DSP_BLOCKMODEL_H
#define DSP_BLOCKMODEL_H

/** standard */
#include <map>
/** Dsp */
#include "Model/DetBlock.h"

using namespace std;

/**
 * A block model is a deterministic model that consists of a block of the original deterministic model. A purpose
 * of this is to represent general decomposable structure, including Benders decomposition and Dantzig-Wolfe
 * decomposition.
 *
 * Each block represents a block of rows. For example, a problem described in DecModel.h has k+1 block models:
 * the first block is
 *   lb^0 <= H^1 x^1 + ... + H^k x^k                     <= ub^0,
 * the second block is
 *   lb^1 <= A^1 x^1                 + B^1 y^1           <= ub^1.
 * and the last block is
 *   lb^k <= A^k x^k                           + B^k y^k <= ub^k.
 *
 * NOTE: These blocks are connected by DecBlkModel.h.
 */

/**
 * Definitions:
 * Primal Block-Angular Matrix: Any two subproblems do not share columns.
 *   X X X X
 *     X
 *       X
 *         X
 * Dual Block-Angular Matrix: Subproblems may share the only columns that are also coupled with the master.
 *   X
 *   X X
 *   X   X
 *   X     X
 */

/**
 * NOTE1: If any two subproblems share columns, the constraint matrix of the full problem is not primal block-angular.
 *        For example,
 *          X X X
 *            X   X
 *            X     X
 *              X     X
 * NOTE2: If the number of columns for every subproblem is less than or equal to that of the master problem,
 *        the constraint matrix of the full problem is not dual block-angular.
 * NOTE3: If any two subproblems share non-coupling columns (i.e., columns not coupled with the master),
 *        the subproblems are not decomposable. But, we assume that they are always decomposable. The column indices
 *        larger than the master's are always decoupled.
 * NOTE4: A constraint matrix may be both primal and dual block-angular. For example,
 *          X   X   X
 *          X X
 *              X X
 *                  X X
 */

typedef map<int,DetBlock*> Blocks;

class BlkModel {

public:

    /** default constructor */
    BlkModel();

    /** copy constructor */
    BlkModel(const BlkModel& rhs);

    /** default destructor */
    virtual ~BlkModel() {}

    /** add block */
    DSP_RTN_CODE addBlock(int id, DetBlock* block);

    /** update coupling columns/rows between the master and the others */
    DSP_RTN_CODE updateBlocks();

    /** get block */
    DetBlock* block(int id);

    /** get blocks */
    Blocks blocks() {return blocks_;}

    /** get number of blocks */
    int getNumBlocks() {return blocks_.size();}

    /** get number of full columns */
    int getNumFullCols() {return ncols_full_;}

    /** get number of full rows */
    int getNumFullRows() {return nrows_full_;}

    /** get number of integers */
    int getNumIntegers() {return nints_full_;}

    vector<int> getCoupledSubproblemIndices(int j) {return coupled_subproblems_indices_[j];}

protected:

    Blocks blocks_; /**< model blocks */

    int nrows_full_; /**< number of rows of the full model */
    int ncols_full_; /**< number of columns of the full model */
    int nints_full_; /**< number of integers of the full model */

    bool primal_block_angular_; /**< Is the full matrix a primal block angular matrix? */
    bool dual_block_angular_;   /**< Is the full matrix a dual block angular matrix? */

    vector<vector<int>> coupled_subproblems_indices_;

};


#endif //DSP_BLOCKMODEL_H
