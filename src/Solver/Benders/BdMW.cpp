/*
 * BdMW.cpp
 *
 *  Created on: Feb 17, 2016
 *      Author: kibaekkim
 */

#include <algorithm>
#include "Solver/Benders/BdMW.h"

BdMW::BdMW(
		DecModel *   model,  /**< model pointer */
		DspParams *  par,    /**< parameters */
		DspMessage * message /**< message pointer */):
BaseMasterWorker(),
model_(model),
par_(par),
message_(message),
master_(NULL),
worker_(NULL) {}

BdMW::BdMW(const BdMW& rhs) :
BaseMasterWorker(rhs),
model_(rhs.model_),
par_(rhs.par_),
message_(rhs.message_) {
	master_ = rhs.master_->clone();
	worker_ = rhs.worker_->clone();
}

BdMW::~BdMW() {
	model_   = NULL;
	par_     = NULL;
	message_ = NULL;
	FREE_PTR(master_);
	FREE_PTR(worker_);
}
