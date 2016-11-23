/*
 * DwMasterTr.h
 *
 *  Created on: Nov 18, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DANTZIGWOLFE_DWMASTERTR_H_
#define SRC_SOLVER_DANTZIGWOLFE_DWMASTERTR_H_

#include <DantzigWolfe/DwMaster.h>

class DwMasterTr: public DwMaster {
public:

    /** default constructor */
	DwMasterTr(
    		DecModel *   model,  /**< model pointer */
            DspParams *  par,    /**< parameters */
            DspMessage * message /**< message pointer */):
	DwMaster(model, par, message),
	tr_cnt_(0),
	tr_size_(0.01),
	tr_center_(NULL) {}

    /** default destructor */
	virtual ~DwMasterTr() {
		FREE_ARRAY_PTR(tr_center_);
	}

    /** solve */
    virtual DSP_RTN_CODE solve();

protected:

	/** Generate initial columns */
	virtual DSP_RTN_CODE initialColumns() {return DSP_RTN_OK;}

    /** This creates a master problem. */
	virtual DSP_RTN_CODE createProblem();

	/** update master */
	virtual DSP_RTN_CODE updateModel(
    		const double* price, /**< [in] price */
			double curLb         /**< [in] current lower bound */);

    /** termination test */
    virtual bool terminationTest(int nnewcols, int itercnt, double relgap);

private:

	/** Check if the current dual is on the trust region boundary */
	bool isTrBoundary(const double* price);

	/** update trust region and model */
	DSP_RTN_CODE updateTrustRegion();

	int tr_cnt_; /**< null step counter */
	double tr_size_; /**< trust region size */
	double* tr_center_; /**< trust region center */
};

#endif /* SRC_SOLVER_DANTZIGWOLFE_DWMASTERTR_H_ */
