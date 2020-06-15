/*
 * DdMaster.cpp
 *
 *  Created on: Feb 9, 2016
 *      Author: kibaekkim
 */

#include "CoinHelperFunctions.hpp"
#include "Model/DecTssModel.h"
#include "Solver/DualDecomp/DdMaster.h"
#include "SolverInterface/DspOsiClp.h"
#include "SolverInterface/DspOsiCpx.h"
#include "SolverInterface/DspOsiGrb.h"
#include "SolverInterface/DspOsiScip.h"
#include "SolverInterface/DspOsiOoqp.h"

DdMaster::DdMaster(DecModel* model, DspParams* par, DspMessage* message) :
		DecSolver(model, par, message),
lambda_(NULL),
subsolution_(NULL) {
	/**< nothing to do */
}

DdMaster::DdMaster(const DdMaster& rhs) :
DecSolver(rhs),
lambda_(rhs.lambda_) {
	subsolution_ = new double * [model_->getNumSubproblems()];
	if (model_->isStochastic()) {
		DecTssModel* tss = dynamic_cast<DecTssModel*>(model_);
		int nsubsolution = tss->getNumCols(0) + tss->getNumCols(1) + 1;
		for (int s = 0; s < model_->getNumSubproblems(); ++s) {
			subsolution_[s] = new double [nsubsolution];
			CoinCopyN(rhs.subsolution_[s], nsubsolution, subsolution_[s]);
		}
	} else {
		for (int s = 0; s < model_->getNumSubproblems(); ++s) {
			subsolution_[s] = new double [model_->getNumSubproblemCouplingCols(s)];
			CoinCopyN(rhs.subsolution_[s], model_->getNumSubproblemCouplingCols(s), subsolution_[s]);
		}
	}
}

DdMaster::~DdMaster() {
	lambda_ = NULL;
	FREE_2D_ARRAY_PTR(model_->getNumSubproblems(),subsolution_);
}

DSP_RTN_CODE DdMaster::init() {
	BGN_TRY_CATCH

	/** status */
	status_ = DSP_STAT_MW_CONTINUE;

	/** time stamp */
	ticToc();

	/** allocate memory */
	subsolution_ = new double * [model_->getNumSubproblems()];
	if (model_->isStochastic()) {
		DecTssModel* tss = dynamic_cast<DecTssModel*>(model_);
		int nsubsolution = tss->getNumCols(0) + tss->getNumCols(1) + 1;
		for (int s = 0; s < model_->getNumSubproblems(); ++s)
			subsolution_[s] = new double [nsubsolution];
	} else {
		for (int s = 0; s < model_->getNumSubproblems(); ++s)
			subsolution_[s] = new double [model_->getNumSubproblemCouplingCols(s)];
	}

	/** initialize */
	primsol_.resize(model_->getNumCouplingCols());
	bestprimsol_.resize(model_->getFullModelNumCols());
	bestdualsol_.resize(model_->getNumCouplingRows());
	subprimobj_.resize(model_->getNumSubproblems());
	subdualobj_.resize(model_->getNumSubproblems());

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}


/** set init solution */
DSP_RTN_CODE DdMaster::setInitSolution(const double * sol) {
	BGN_TRY_CATCH

	if (primsol_.size() >= getSiPtr()->getNumCols())
		CoinCopyN(sol, getSiPtr()->getNumCols(), &primsol_[0]);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DspOsi * DdMaster::createDspOsi() {
	DspOsi * osi = NULL;
	BGN_TRY_CATCH

	switch (par_->getIntParam("DD/MASTER/SOLVER")) {
	case OsiCpx:
#ifdef DSP_HAS_CPX
		osi = new DspOsiCpx();
#else
		throw CoinError("Cplex is not available.", "createDspOsi", "DdMaster");
#endif
		break;
	case OsiGrb:
#ifdef DSP_HAS_GRB
		osi = new DspOsiGrb();
#else
		throw CoinError("Gurobi is not available.", "createDspOsi", "DdMaster");
#endif
		break;
	case OsiScip:
#ifdef DSP_HAS_SCIP
		osi = new DspOsiScip();
#else
		throw CoinError("Scip is not available.", "createDspOsi", "DdMaster");
#endif
		break;
	case OsiOoqp:
#ifdef DSP_HAS_OOQP
		osi = new DspOsiOoqp();
#else
		throw CoinError("OOQP is not available.", "createDspOsi", "DdMaster");
#endif
		break;
	case OsiClp:
		osi = new DspOsiClp();
		break;
	default:
		throw CoinError("Invalid parameter value", "createDspOsi", "DdMaster");
		break;
	}

	END_TRY_CATCH(;)
	
	return osi;
}
