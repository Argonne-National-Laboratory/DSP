/*
 * DwBundleDualPips.h
 *
 *  Created on: Dec 7, 2017
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DANTZIGWOLFE_DWBUNDLEDUALPIPS_H_
#define SRC_SOLVER_DANTZIGWOLFE_DWBUNDLEDUALPIPS_H_

#include "Solver/DantzigWolfe/DwBundleDual.h"

class DwBundleDualPips: public DwBundleDual {
public:
	/** constructor with worker */
	DwBundleDualPips(DwWorker* worker) : DwBundleDual(worker) {}

	/** copy constructor */
	DwBundleDualPips(const DwBundleDualPips& rhs) : DwBundleDual(rhs) {
	}

	/** copy operator */
	DwBundleDualPips& operator=(const DwBundleDualPips& rhs);

	virtual DwBundleDualPips* clone() const {
		return new DwBundleDualPips(*this);
	}

	/** default destructor */
	virtual ~DwBundleDualPips() {}

protected:

	/** Update the proximal bundle center */
	virtual DSP_RTN_CODE updateCenter(double penalty);

	/** print iteration information */
	virtual void printIterInfo();

	//@{
	/** functions specific to external solver */

	/** initialize dual solver */
	virtual void initDualSolver(
			const CoinPackedMatrix& m, 
			std::vector<double>& clbd, 
			std::vector<double>& cubd, 
			std::vector<double>& obj);

	/** call master solver */
	virtual DSP_RTN_CODE callMasterSolver();

	/** get master solution */
	virtual void assignMasterSolution(std::vector<double>& sol);

	/** get objective value */
	virtual double getObjValue();

	/** add row to the dual problem */
	virtual void addDualRow(const CoinPackedVector& v, const double lb, const double ub);

	/** remove all rows in the dual master */
	virtual void removeAllDualRows();

	//@}
};

#endif /* SRC_SOLVER_DANTZIGWOLFE_DWBUNDLEDUAL_H_ */
