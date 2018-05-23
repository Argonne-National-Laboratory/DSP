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
		active_(true),
		age_(0) {
		/** that's it */
	}

	/** constructor */
	DwCol(int blockid, CoinPackedVector x, CoinPackedVector col, double obj, double lb, double ub, bool active):
		blockid_(blockid),
		x_(x),
		col_(col),
		obj_(obj),
		lb_(lb),
		ub_(ub),
		active_(active),
		age_(0) {
		/** that's it */
	}

	/** copy constructor */
	DwCol(const DwCol& rhs):
		blockid_(rhs.blockid_),
		x_(rhs.x_),
		col_(rhs.col_),
		obj_(rhs.obj_),
		lb_(rhs.lb_),
		ub_(rhs.ub_),
		active_(rhs.active_),
		age_(rhs.age_) {
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
	int age_;
};



#endif /* SRC_SOLVER_DANTZIGWOLFE_DWCOL_H_ */
