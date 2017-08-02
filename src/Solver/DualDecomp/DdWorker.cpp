/*
 * DdWorker.cpp
 *
 *  Created on: Feb 9, 2016
 *      Author: kibaekkim
 */

/** DSP */
#include "Solver/DualDecomp/DdWorker.h"

DdWorker::DdWorker(DspParams * par, DecModel * model, DspMessage * message) :
		DdSolver(par, model, message) {
}

DdWorker::~DdWorker() {
}
