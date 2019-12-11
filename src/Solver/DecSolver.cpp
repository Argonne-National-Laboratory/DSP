/*
 * DecSolver.cpp
 *
 *  Created on: Sep 24, 2014
 *      Author: kibaekkim, ctjandra
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include "Solver/DecSolver.h"

#ifdef DSP_HAS_CPX
#include "OsiCpxSolverInterface.hpp"
#endif

#ifdef DSP_HAS_SCIP
#include "SolverInterface/OsiScipSolverInterface.hpp"
#endif

DecSolver::DecSolver(
			DecModel *   model,   /**< model pointer */
			DspParams *  par,     /**< parameter pointer */
			DspMessage * message /**< message pointer */) :
model_(model),
par_(par),
message_(message),
si_(NULL),
status_(DSP_STAT_UNKNOWN),
bestprimobj_(COIN_DBL_MAX),
primobj_(COIN_DBL_MAX),
bestdualobj_(-COIN_DBL_MAX),
dualobj_(-COIN_DBL_MAX),
cputime_(0.0),
walltime_(0.0),
time_remains_(COIN_DBL_MAX),
tic_(0.0),
numIterations_(0),
numNodes_(0),
iterlim_(COIN_INT_MAX) {}

DecSolver::DecSolver(const DecSolver&rhs) :
model_(rhs.model_),
par_(rhs.par_),
message_(rhs.message_),
status_(rhs.status_),
bestprimobj_(rhs.bestprimobj_),
primobj_(rhs.primobj_),
bestdualobj_(rhs.bestdualobj_),
dualobj_(rhs.dualobj_),
cputime_(rhs.cputime_),
walltime_(rhs.walltime_),
time_remains_(rhs.time_remains_),
tic_(rhs.tic_),
numIterations_(rhs.numIterations_),
numNodes_(rhs.numNodes_),
iterlim_(rhs.iterlim_) {
	si_ = rhs.si_->clone();
	bestprimsol_ = rhs.bestprimsol_;
	primsol_ = rhs.primsol_;
	bestdualsol_ = rhs.bestdualsol_;
	dualsol_ = rhs.dualsol_;
	s_statuses_ = rhs.s_statuses_;
	s_primobjs_ = rhs.s_primobjs_;
	s_dualobjs_ = rhs.s_dualobjs_;
	s_primsols_ = rhs.s_primsols_;
	s_cputimes_ = rhs.s_cputimes_;
	s_walltimes_ = rhs.s_walltimes_;
}

DecSolver::~DecSolver() {
	FREE_PTR(si_);
	message_ = NULL;
	par_ = NULL;
	model_ = NULL;
}

/** solution status */
DSP_RTN_CODE DecSolver::getStatus()
{
	return getStatus(si_);
}

/** solution status */
DSP_RTN_CODE DecSolver::getStatus(OsiSolverInterface* si)
{
	int status = DSP_STAT_UNKNOWN;
	if (si->isProvenOptimal())
		status = DSP_STAT_OPTIMAL;
	else if (si->isProvenDualInfeasible())
		status = DSP_STAT_DUAL_INFEASIBLE;
	else if (si->isProvenPrimalInfeasible())
		status = DSP_STAT_PRIM_INFEASIBLE;
	else if (si->isPrimalObjectiveLimitReached())
		status = DSP_STAT_LIM_PRIM_OBJ;
	else if (si->isDualObjectiveLimitReached())
		status = DSP_STAT_LIM_DUAL_OBJ;
	else if (si->isIterationLimitReached())
		status = DSP_STAT_LIM_ITERorTIME;
	else if (si->isAbandoned())
		status = DSP_STAT_STOPPED_UNKNOWN;
	return status;
}

/** write output to a file */
void DecSolver::write(const char * filename)
{
	BGN_TRY_CATCH

	ofstream myfile;
	myfile.open(filename);
	myfile << "Iter";
	myfile << ",Status";
	myfile << ",Prim";
	myfile << ",Dual";
	myfile << ",Cpu";
	myfile << ",Wall";
	myfile << "\n";
	for (unsigned i = 0; i < s_statuses_.size(); ++i)
	{
		myfile << i;
		myfile << "," << s_statuses_[i];
		myfile << "," << scientific << setprecision(5) << s_primobjs_[i];
		myfile << "," << scientific << setprecision(5) << s_dualobjs_[i];
		myfile << "," << fixed << setprecision(2) << s_cputimes_[i];
		myfile << "," << fixed << setprecision(2) << s_walltimes_[i];
		myfile << "\n";
	}
	myfile.close();

	END_TRY_CATCH(;)
}
