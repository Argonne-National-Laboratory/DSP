/*
 * DwBundleDualSmip.h
 *
 *  Created on: Jun 27, 2018
 *      Author: Kibaek Kim
 */

#ifndef SRC_SOLVER_DANTZIGWOLFE_DWBUNDLEDUALSMIP_H_
#define SRC_SOLVER_DANTZIGWOLFE_DWBUNDLEDUALSMIP_H_

#include <Model/TssModel.h>
#include <DantzigWolfe/DwBundleDual.h>

class DwBundleDualSmip: public DwBundleDual {
public:
	/** constructor with worker */
	DwBundleDualSmip(DwWorker* worker);

	/** copy constructor */
	DwBundleDualSmip(const DwBundleDualSmip& rhs);

	/** copy operator */
	DwBundleDualSmip& operator=(const DwBundleDualSmip& rhs);

	virtual DwBundleDualSmip* clone() const {
		return new DwBundleDualSmip(*this);
	}

	/** default destructor */
	virtual ~DwBundleDualSmip();

    /** initialize */
    virtual DSP_RTN_CODE init();

protected:

	/** create primal master problem */
	virtual DSP_RTN_CODE createPrimalProblem();

	/** create dual master problem */
	virtual DSP_RTN_CODE createDualProblem();

	//@{
	/** functions specific to external solver */

	/** remove all columns in the primal master */
	virtual void removeAllPrimCols();

	/** remove all rows in the dual master */
	virtual void removeAllDualRows();

	//@}

	TssModel* tss_;
};

#endif /* SRC_SOLVER_DANTZIGWOLFE_DWBUNDLEDUALSMIP_H_ */
