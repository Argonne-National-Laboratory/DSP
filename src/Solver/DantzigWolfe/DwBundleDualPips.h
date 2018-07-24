/*
 * DwBundleDualPips.h
 *
 *  Created on: Dec 7, 2017
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DANTZIGWOLFE_DWBUNDLEDUALPIPS_H_
#define SRC_SOLVER_DANTZIGWOLFE_DWBUNDLEDUALPIPS_H_

#include "Solver/DantzigWolfe/DwWorkerPips.h"
#include "Solver/DantzigWolfe/DwBundleDualSmip.h"

class DwBundleDualPips: public DwBundleDualSmip {
public:
	/** constructor with worker */
	DwBundleDualPips(DwWorker* worker) : DwBundleDualSmip(worker) {
		pips_worker_ = dynamic_cast<DwWorkerPips*>(worker_);
	}

	/** copy constructor */
	DwBundleDualPips(const DwBundleDualPips& rhs) : DwBundleDualSmip(rhs) {
		pips_worker_ = rhs.pips_worker_;
	}

	/** copy operator */
	DwBundleDualPips& operator=(const DwBundleDualPips& rhs);

	virtual DwBundleDualPips* clone() const {
		return new DwBundleDualPips(*this);
	}

	/** default destructor */
	virtual ~DwBundleDualPips() {
		pips_worker_ = NULL;
	}

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
			std::vector<double>& obj, 
			std::vector<double>& rlbd, 
			std::vector<double>& rubd);

	/** call master solver */
	virtual DSP_RTN_CODE callMasterSolver();

	/** get master solution */
	virtual void assignMasterSolution(std::vector<double>& sol);

	/** get objective value */
	virtual double getObjValue();

	/** add row to the dual problem */
	virtual void addDualRow(const CoinPackedVector& v, double lb, double ub);

	/** remove all rows in the dual master */
	virtual void removeAllDualRows();

	//@}

	/** pips worker */
	DwWorkerPips* pips_worker_;
};

#endif /* SRC_SOLVER_DANTZIGWOLFE_DWBUNDLEDUAL_H_ */
