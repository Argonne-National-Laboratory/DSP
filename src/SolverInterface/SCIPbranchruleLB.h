/*
 * SCIPbranchruleLB.h
 *
 *  Created on: Apr 9, 2015
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_SCIPBRANCHRULELB_H_
#define SRC_SOLVER_SCIPBRANCHRULELB_H_

/** SCIP */
#include "objscip/objbranchrule.h"

#include "CoinFinite.hpp"

/** This is to prune nodes based on a known global lower bound */
class SCIPbranchruleLB : public scip::ObjBranchrule {
public:

	/** default constructor */
	SCIPbranchruleLB(SCIP * scip) :
		ObjBranchrule(scip, "knownLB", "knownLB", 100000, -1, 1.0), lowerbound_(-COIN_DBL_MAX)
	{
		/** TODO */
	}

	virtual ~SCIPbranchruleLB()
	{
		/** TODO */
	}

	/** branching execution method for fractional LP solutions
	 *
	 *  @see SCIP_DECL_BRANCHEXECLP(x) in @ref type_branch.h
	 */
	virtual SCIP_DECL_BRANCHEXECLP(scip_execlp);

	/** branching execution method for external candidates
	 *
	 *  @see SCIP_DECL_BRANCHEXECEXT(x) in @ref type_branch.h
	 */
	virtual SCIP_DECL_BRANCHEXECEXT(scip_execext);

	/** branching execution method for not completely fixed pseudo solutions
	 *
	 *  @see SCIP_DECL_BRANCHEXECPS(x) in @ref type_branch.h
	 */
	virtual SCIP_DECL_BRANCHEXECPS(scip_execps);

	/** set global lower bound */
	void setLowerBound(double lb) {lowerbound_ = lb;}

	/** check lower bound */
	void checkLowerBound(SCIP * scip, SCIP_RESULT * result);

private:

	double lowerbound_;
};

#endif /* SRC_SOLVER_SCIPBRANCHRULELB_H_ */
