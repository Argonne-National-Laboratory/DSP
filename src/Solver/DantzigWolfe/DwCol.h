/*
 * DwCol.h
 *
 *  Created on: Oct 11, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DANTZIGWOLFE_DWCOL_H_
#define SRC_SOLVER_DANTZIGWOLFE_DWCOL_H_

class DwCol {
public:

	/** constructor */
	DwCol(int blockid, CoinPackedVector x, CoinPackedVector col, double obj, double lb, double ub):
		blockid_(blockid),
		x_(x),
		col_(col),
		obj_(obj),
		lb_(lb),
		ub_(ub),
		active_(true) {
		/** that's it */
	}

	/** destructor */
	virtual ~DwCol() {}

	int blockid_; /**< subproblem block id based on zero */
	CoinPackedVector x_;
	CoinPackedVector col_;
	double obj_;
	double lb_;
	double ub_;
	bool active_;
};



#endif /* SRC_SOLVER_DANTZIGWOLFE_DWCOL_H_ */
