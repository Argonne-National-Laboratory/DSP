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
