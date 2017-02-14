/*
 * DwModel.h
 *
 *  Created on: Feb 13, 2017
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DANTZIGWOLFE_DWMODEL_H_
#define SRC_SOLVER_DANTZIGWOLFE_DWMODEL_H_

#include <TreeSearch/DspModel.h>

class DwModel: public DspModel {
public:
	/** default constructor */
	DwModel();

	/** default constructor with solver */
	DwModel(DecSolver* solver);

	/** default destructor */
	virtual ~DwModel();

	/** solve model */
    virtual DSP_RTN_CODE solve();
};

#endif /* SRC_SOLVER_DANTZIGWOLFE_DWMODEL_H_ */
