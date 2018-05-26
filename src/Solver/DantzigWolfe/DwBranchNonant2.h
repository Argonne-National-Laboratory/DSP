/*
 * DwBranchNonant2.h
 *
 *  Created on: May 8, 2018
 *      Author: Kibaek Kim
 * 
 *  This implements the branching on the nonanticipativity constraints, generalized from that of Caroe and Schultz (1999).
 */

#ifndef SRC_SOLVER_DANTZIGWOLFE_DWBRANCHNONANT2_H_
#define SRC_SOLVER_DANTZIGWOLFE_DWBRANCHNONANT2_H_

#include <DantzigWolfe/DwBranchNonant.h>

class DwBranchNonant2 : public DwBranchNonant {
public:
	/** default constructor */
	DwBranchNonant2() : DwBranchNonant() {}

	/** default constructor with solver */
	DwBranchNonant2(DwModel* model) : DwBranchNonant(model) {}

	/** default destructor */
	virtual ~DwBranchNonant2() {
	}

    /** calculate and return reference solution */
    virtual void getRefSol(std::vector<double>& refsol);

    /** calculate and return deviation of the reference solution */
    virtual void getDevSol(std::vector<double>& refsol, std::vector<double>& devsol);

};

#endif /* SRC_SOLVER_DANTZIGWOLFE_DWBRANCHNONANT_H_ */
