/*
 * BdDriver.cpp
 *
 *  Created on: Feb 17, 2016
 *      Author: kibaekkim
 */

#include "Solver/Benders/BdDriver.h"

BdDriver::BdDriver(
		DecModel *   model,  /**< model pointer */
		DspParams *  par,    /**< parameters */
		DspMessage * message /**< message pointer */):
DecSolver(model, par, message),
mw_(NULL),
aux_size_(0),
aux_obj_(NULL),
aux_clbd_(NULL),
aux_cubd_(NULL),
numPriorities_(0),
priorities_(NULL) {
	par_->setDblParam("DD/WALL_LIM", par_->getDblParam("BD/WALL_LIM"));
}

BdDriver::BdDriver(const BdDriver&rhs) :
DecSolver(rhs),
aux_size_(rhs.aux_size_),
numPriorities_(rhs.numPriorities_) {
	// copy master-worker framework
	mw_ = rhs.mw_->clone();

	// copy auxiliary variable
	setAuxVarData(rhs.aux_size_, rhs.aux_obj_, rhs.aux_clbd_, rhs.aux_cubd_);

	// copy priority
	priorities_ = NULL;
	setPriorities(rhs.numPriorities_, rhs.priorities_);
}

BdDriver::~BdDriver() {
	FREE_PTR(mw_);
	FREE_ARRAY_PTR(aux_obj_);
	FREE_ARRAY_PTR(aux_clbd_);
	FREE_ARRAY_PTR(aux_cubd_);
	FREE_ARRAY_PTR(priorities_);
}

DSP_RTN_CODE BdDriver::setAuxVarData(
		int      size, /**< size of arrays */
		double * obj,  /**< objective function coefficients */
		double * clbd, /**< column lower bounds */
		double * cubd  /**< column upper bounds */)
{
	if (size <= 0) return DSP_RTN_ERR;

	BGN_TRY_CATCH

	FREE_ARRAY_PTR(aux_obj_)
	FREE_ARRAY_PTR(aux_clbd_)
	FREE_ARRAY_PTR(aux_cubd_)

	aux_size_ = size;
	aux_obj_ = new double [size];
	aux_clbd_ = new double [size];
	aux_cubd_ = new double [size];
	CoinCopyN(obj, size, aux_obj_);
	CoinCopyN(clbd, size, aux_clbd_);
	CoinCopyN(cubd, size, aux_cubd_);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE BdDriver::setPriorities(
		int   size,      /**< size of array */
		int * priorities /**< branch priority */)
{
	if (size <= 0) return DSP_RTN_ERR;

	BGN_TRY_CATCH

	numPriorities_ = size;
	if (priorities_ == NULL)
		priorities_ = new int [numPriorities_];
	CoinCopyN(priorities, numPriorities_, priorities_);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE BdDriver::setSolution(
		int      size,    /**< size of array */
		double * solution /**< solution */)
{
	if (size <= 0) return DSP_RTN_ERR;

	BGN_TRY_CATCH

	initsols_.push_back(new CoinPackedVector(size, solution));

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}
