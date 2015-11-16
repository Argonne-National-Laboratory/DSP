/*
 * DecTssSolver.h
 *
 *  Created on: Jul 2, 2015
 *      Author: ctjandra
 */

#ifndef DECTSSSOLVER_H_
#define DECTSSSOLVER_H_

/** DSP */
#include "Utility/StoConfig.h"
#include "Utility/StoRtnCodes.h"
#include "Utility/StoUtility.h"
#include "Utility/StoMessage.h"
#include "Utility/DspParams.h"
#include "Model/DecModel.h"

/** 
 * TODO:
 * This is a TEMPORARY class that disguises a TssSolver as a DecSolver.
 * This class should be deleted once all TssSolvers are migrated to DecSolvers.
 */

class DecTssSolver : public DecSolver {
public:

	/** default constructor */
	DecTssSolver() : DecSolver(), tss_(NULL) {}

	DecTssSolver(TssSolver * tss) : tss_(tss) {}

	/** default destructor */
	virtual ~DecTssSolver() {
		solution_ = NULL;
		FREE_PTR(tss_);
	}

	/** load model object; will have a shallow pointer (not deep copy) */
	STO_RTN_CODE loadModel(DspParams * par, DecModel * model)
	{
		TssModel * tssModel = dynamic_cast<TssModel *>(model);
		if (tssModel == NULL)
		{
			printf("Error: Model in DecTssSolver must be TssModel");
			return STO_RTN_ERR;
		}
		return tss_->loadModel(par, tssModel);
	}

	/** solve */
	STO_RTN_CODE solve()
	{
		int rtn_code = tss_->solve();
		status_ = tss_->status_;
		solution_ = tss_->solution_;
		primalBound_ = tss_->primalBound_;
		dualBound_ = tss_->dualBound_;
		numIterations_ = tss_->numIterations_;
		numNodes_ = tss_->numNodes_;
		clockType_ = tss_->clockType_;
		solutionTime_ = tss_->solutionTime_;
		return rtn_code;
	}

public:

	TssSolver * tss_;
};

#endif /* DECTSSSOLVER_H_ */
