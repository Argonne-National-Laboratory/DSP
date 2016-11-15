/*
 * BdMaster.cpp
 *
 *  Created on: Feb 18, 2016
 *      Author: kibaekkim
 */

//#define DSP_DEBUG

/** Coin */
#include "OsiClpSolverInterface.hpp"
/** Dsp */
#include "Solver/Benders/BdMaster.h"

BdMaster::BdMaster(DecModel* model, DspParams* par, DspMessage * message) :
		DecSolver(model, par, message),
		si_(NULL),
		naux_(1),
		obj_aux_(NULL),
		clbd_aux_(NULL),
		cubd_aux_(NULL),
		worker_(NULL) {
	obj_aux_ = new double[naux_];
	clbd_aux_ = new double[naux_];
	cubd_aux_ = new double[naux_];
	obj_aux_[0] = 1.0;
	clbd_aux_[0] = -COIN_DBL_MAX;
	cubd_aux_[0] = +COIN_DBL_MAX;
	tic_ = CoinGetTimeOfDay();
}

BdMaster::~BdMaster() {
	FREE_PTR(si_);
	FREE_ARRAY_PTR(obj_aux_);
	FREE_ARRAY_PTR(clbd_aux_);
	FREE_ARRAY_PTR(cubd_aux_);
}

DSP_RTN_CODE BdMaster::init() {
	BGN_TRY_CATCH

	/** create problem */
	createProblem();

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE BdMaster::solve() {

	/** warm start information */
	CoinWarmStartBasis* ws = NULL;

	BGN_TRY_CATCH

	DSPdebugMessage("Start solving...\n");

	/** initial solve */
	si_->resolve();

	/** status */
	convertCoinToDspStatus(si_, status_);

	OsiCuts cuts;
	while (status_ == DSP_STAT_OPTIMAL) {
		/** generate cuts */
		worker_->generateCuts(si_->getNumCols() - naux_, naux_, si_->getColSolution(), cuts);

		/** evaluate cuts */
		int nCutsToAdd = 0;
		double infeasibility = 0.0;
		for (int i = 0; i < cuts.sizeRowCuts(); ++i) {
			double eff = cuts.rowCutPtr(i)->violated(si_->getColSolution());
			cuts.rowCutPtr(i)->setEffectiveness(eff);
			infeasibility += eff;
			if (eff > 1.0e-6) nCutsToAdd++;
		}

		/** done? */
		if (nCutsToAdd == 0) break;

		/** get basis */
		FREE_PTR(ws);
		ws = dynamic_cast<CoinWarmStartBasis*>(si_->getWarmStart());
		int nrows = si_->getNumRows();
		ws->resize(nrows + nCutsToAdd, si_->getNumCols());

		/** add cuts */
		for (int i = 0; i < cuts.sizeRowCuts(); ++i) {
			if (cuts.rowCutPtr(i)->effectiveness() > 1.0e-6)
				si_->applyRowCut(cuts.rowCut(i));
		}

		/** adjust basis */
		for (int i = 0; i < nCutsToAdd; ++i)
			ws->setArtifStatus(nrows + i, CoinWarmStartBasis::basic);
		si_->setWarmStart(ws);

		/** re-optimize the master */
		si_->resolve();
		convertCoinToDspStatus(si_, status_);
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE BdMaster::setSolutions(Solutions initsols) {
	throw CoinError("Not implemented.", "setSolutions", "BdMaster");
	return DSP_RTN_OK;
}

DSP_RTN_CODE BdMaster::finalize() {
	return DSP_RTN_OK;
}

DSP_RTN_CODE BdMaster::createProblem() {
#define FREE_MEMORY       \
	FREE_PTR(mat)         \
	FREE_ARRAY_PTR(clbd)  \
	FREE_ARRAY_PTR(cubd)  \
	FREE_ARRAY_PTR(ctype) \
	FREE_ARRAY_PTR(obj)   \
	FREE_ARRAY_PTR(rlbd)  \
	FREE_ARRAY_PTR(rubd)

	assert(model_);

	if (naux_ <= 0 || !obj_aux_ || !clbd_aux_ || !cubd_aux_) {
		printf("Warning: Auxiliary column information is required.\n");
		return DSP_RTN_ERR;
	}

	/** master problem */
	CoinPackedMatrix * mat = NULL;
	double * clbd   = NULL;
	double * cubd   = NULL;
	double * obj    = NULL;
	char *   ctype  = NULL;
	double * rlbd   = NULL;
	double * rubd   = NULL;

	BGN_TRY_CATCH

	/** number of columns */
	int ncols = model_->getNumCouplingCols() + naux_;

	/** decompose model */
	DSP_RTN_CHECK_THROW(model_->decompose(
			0, NULL, naux_, clbd_aux_, cubd_aux_, obj_aux_,
			mat, clbd, cubd, ctype, obj, rlbd, rubd));
	DSPdebugMessage("Decomposed the model.\n");

	assert(si_==NULL);
	si_ = new OsiClpSolverInterface();
	si_->messageHandler()->logLevel(par_->getIntParam("LOG_LEVEL"));
	DSPdebugMessage("Successfully created SCIP interface \n");

	/** load problem data */
	si_->loadProblem(*mat, clbd, cubd, obj, rlbd, rubd);
	DSPdebugMessage("Successfully load the problem.\n");

	/** set column type */
	for (int j = 0; j < ncols; ++j) {
		if (ctype[j] != 'C')
			si_->setInteger(j);
	}

	/** allocate memory for primal solution */
	primsol_ = new double [si_->getNumCols()];

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	/** save memory */
	FREE_MEMORY

	return DSP_RTN_OK;

#undef FREE_MEMORY
}

DSP_RTN_CODE BdMaster::setObjectiveBounds(double upper, double lower) {
#define FREE_MEMORY         \
	FREE_ARRAY_PTR(auxind)  \
	FREE_ARRAY_PTR(auxcoef)

	/** for recourse lower bound */
	int * auxind     = NULL;
	double * auxcoef = NULL;

	BGN_TRY_CATCH

	/** number of columns */
	int ncols = si_->getNumCols();

	/** allocate memory */
	auxind = new int [ncols];
	auxcoef = new double [ncols];

	/** update bounds */
	primobj_ = upper;
	dualobj_ = lower;

	for (int j = 0; j < ncols; ++j) {
		auxind[j] = j;
		auxcoef[j] = si_->getObjCoefficients()[j];
	}
	si_->addRow(ncols, auxind, auxcoef, dualobj_, primobj_);

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

DSP_RTN_CODE BdMaster::setAuxVarData(
		int size,
		double* obj,
		double* clbd,
		double* cubd) {
	BGN_TRY_CATCH

	/** free memory (just in case) */
	FREE_ARRAY_PTR(obj_aux_)
	FREE_ARRAY_PTR(clbd_aux_)
	FREE_ARRAY_PTR(cubd_aux_)
	/** allocate memory */
	naux_ = size;
	obj_aux_ = new double[naux_];
	clbd_aux_ = new double[naux_];
	cubd_aux_ = new double[naux_];
	/** copy data */
	CoinCopyN(obj, naux_, obj_aux_);
	CoinCopyN(clbd, naux_, clbd_aux_);
	CoinCopyN(cubd, naux_, cubd_aux_);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}
