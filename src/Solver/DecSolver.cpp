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

DecSolver::DecSolver(DspParams * par, DecModel * model, DspMessage * message):
	model_(model), par_(par), message_(message),
	status_(DSP_STAT_UNKNOWN), primsol_(NULL), dualsol_(NULL),
	primobj_(COIN_DBL_MAX), dualobj_(-COIN_DBL_MAX),
	time_remains_(COIN_DBL_MAX), tic_(0.0) {}

DecSolver::~DecSolver()
{
	FREE_ARRAY_PTR(primsol_);
	for (unsigned i = 0; i < s_primsols_.size(); ++i)
		FREE_ARRAY_PTR(s_primsols_[i]);
	message_ = NULL;
	par_ = NULL;
	model_ = NULL;
}

/** update time stamp and time remains */
DSP_RTN_CODE DecSolver::ticToc()
{
	BGN_TRY_CATCH

	time_remains_ -= CoinGetTimeOfDay() - tic_;
	tic_ = CoinGetTimeOfDay();

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

/** write output to a file */
void DecSolver::write(const char * filename)
{
	BGN_TRY_CATCH

	ofstream myfile;
	myfile.open(filename);
	myfile << "Priaml Objective: " << primobj_ << "\n";
	myfile << "Dual Objective: " << dualobj_ << "\n";
	myfile << setw(6) << right << "Iter";
	myfile << setw(10) << "Status";
	myfile << setw(14) << "Prim";
	myfile << setw(14) << "Dual";
	myfile << setw(10) << "Cpu";
	myfile << setw(10) << "Wall";
	myfile << "\n";
	for (unsigned i = 0; i < s_statuses_.size(); ++i)
	{
		myfile << setw(6) << i;
		myfile << setw(10) << s_statuses_[i];
		myfile << setw(14) << scientific << setprecision(5) << s_primobjs_[i];
		myfile << setw(14) << scientific << setprecision(5) << s_dualobjs_[i];
		myfile << setw(10) << fixed << setprecision(2) << s_cputimes_[i];
		myfile << setw(10) << fixed << setprecision(2) << s_walltimes_[i];
		myfile << "\n";
	}
	myfile.close();

	END_TRY_CATCH(;)
}
