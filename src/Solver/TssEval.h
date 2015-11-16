/*
 * TssEval.h
 *
 *  Created on: Mar 16, 2015
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_TSSEVAL_H_
#define SRC_SOLVER_TSSEVAL_H_

#include "Solver/TssSolver.h"

class TssEval: public TssSolver
{
public:

	/** default constructor */
	TssEval();

	/** default destructor */
	virtual ~TssEval();

	/** solve */
	virtual STO_RTN_CODE solve();

	/** set solution to evaluate */
	void setSolution(double * solution);

private:

	bool hasSolution_; /** indicating whether a solution is set or not */

protected:

	/** parameters */
	int parNumCores_;
};

#endif /* SRC_SOLVER_TSSEVAL_H_ */
