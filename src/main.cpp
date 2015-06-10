/*
 * main.cpp
 *
 *  Created on: Nov 6, 2014
 *      Author: Kibaek Kim
 */

#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <limits>

#define USE_MPI
#ifdef USE_MPI
/** MPI */
#include "mpi.h"
#endif

/** COIN packages */
#include "OsiClpSolverInterface.hpp"

/** DSP */
#include "Utility/StoMacros.h"
#include "Utility/StoRtnCodes.h"
#include "Model/TssModel.h"
#include "Solver/StoParam.h"
#include "Solver/TssDe.h"
#include "Solver/TssBd.h"
//#include "Solver/TssDd.h"
#ifdef USE_MPI
#include "Solver/TssDdMpi.h"
#endif

using namespace std;

/** This is an executable DSP from reading SMPS file.
 */

int main(int argc, char ** argv)
{
	if (argc == 1)
	{
		cout << "Usage: dsp -f <smps> [options]\n"
				<< "  Options:\n"
				<< "    -a  algorithm (default: 1)\n"
				<< "        0: Deterministic equivalent (DE)\n"
				<< "        1: Benders decomposition (BD)\n"
				<< "        2: Dual decomposition (DD)\n"
				<< "    -n  [BD] number of cores used (default: 1)\n"
				<< "    -m  [DD] master problem solver type (default: 0 [Clp])\n"
				<< "        0: Clp, 1: OOQP, 2: OOQP-Suboptimal\n"
//				<< "    -t  [DD] enables trust region (default: 1)\n"
//				<< "        0: no, 1: yes\n"
				<< "    -b  [DD] enables feasibility recovery (default: 0)\n"
				<< "        0: no, 1: yes\n"
				<< "    -u  [DD] frequency of feasibility recovery (default: -1)\n"
				<< "        < 0: do not recover feasible solution\n"
				<< "          0: only at the first iteration\n"
				<< "        < 0: frequency of recovery\n"
				<< "    -W  wallclock time limit\n"
				<< "    -N  Tree search node limit (default: no limit)\n"
				<< "    -I  Iteration limit (default: no limit)\n"
				<< "    -p  print level (default: 0)\n"
				<< "    -L  prefix of logfile name\n"
				<< endl;
		return 0;
	}

	char smpsfile[1024];
	int algo       = 1;
	int ncores     = 1;
	int loglevel   = 0;
	int feasibleRecovery = -1;
	int use_ubcuts = 0;
	int useTrustRegion = 1;
	int solverType = 0; /**< Lagrangian master solver type */
	double wall_limit = numeric_limits<double>::max();
	int iter_limit    = numeric_limits<int>::max();
	int node_limit    = -1;
	int benders_aggressive = 1;
	string logfile_prefix;

	int c;
	while ((c = getopt(argc, argv, "f:a:m:n:p:u:b:t:L:W:N:I:")) != -1)
	{
		switch (c)
		{
		case 'f':
			strcpy(smpsfile, optarg);
			break;
		case 'a':
			algo = atoi(optarg);
			break;
		case 'm':
			solverType = atoi(optarg);
			break;
		case 'n':
			ncores = atoi(optarg);
			break;
		case 'p':
			loglevel = atoi(optarg);
			break;
		case 'u':
			feasibleRecovery = atoi(optarg);
			break;
		case 'b':
			use_ubcuts = atoi(optarg);
			break;
		case 't':
			useTrustRegion = atoi(optarg);
			break;
		case 'L':
			logfile_prefix.assign(optarg);
			break;
		case 'W':
			wall_limit = atof(optarg);
			break;
		case 'N':
			node_limit = atoi(optarg);
			break;
		case 'I':
			iter_limit = atoi(optarg);
			break;
		case '?':
			break;
		default:
			break;
		}
	}

#ifdef USE_MPI
	if (algo == 2) MPI_Init(&argc, &argv);
#endif

	TssModel * model = NULL;
	StoParam * par = NULL;
	TssSolver * solver = NULL;

	model = new TssModel;

	int rtn = model->readSmps(smpsfile);
	if (rtn == STO_RTN_ERR)
	{
		if (model)
		{
			delete model;
			model = NULL;
		}
#ifdef USE_MPI
		if (algo == 2) MPI_Finalize();
#endif
		return 0;
	}

	par = new StoParam;
	par->TssBdBendersPriority_ = benders_aggressive;
	par->TssDdEnableTrustRegion_   = useTrustRegion;
	par->TssDdTrustRegionSize_ = 1.;
	par->TssDdMasterSolver_ = solverType;
	par->wtimeLimit_ = wall_limit;
	par->iterLimit_ = iter_limit;
	par->logLevel_ = loglevel;
	par->nodeLimit_ = node_limit;
	par->numCores_ = ncores;

	switch(algo)
	{
	case 0: /** De */
		solver = new TssDe;
		break;
	case 1:
	{
		solver = new TssBd;
		/** set auxiliary variables */
		double * obj_aux = new double [1];
		double * clbd_aux = new double [1];
		double * cubd_aux = new double [1];
		CoinFillN(obj_aux, 1, 1.0);
		CoinFillN(clbd_aux, 1, -COIN_DBL_MAX);
		CoinFillN(cubd_aux, 1, +COIN_DBL_MAX);
		dynamic_cast<TssBd*>(solver)->setAuxColData(1, obj_aux, clbd_aux, cubd_aux);
		delete [] obj_aux;
		delete [] clbd_aux;
		delete [] cubd_aux;
		obj_aux = NULL;
		clbd_aux = NULL;
		cubd_aux = NULL;
		break;
	}
	case 2:
	{
#ifdef USE_MPI
		solver = new TssDdMpi(MPI_COMM_WORLD, logfile_prefix);
#else
		solver = new TssDd;
#endif
		break;
	}
	default:
		break;
	}

	//solver->messageHandler()->setLogLevel(loglevel);
	solver->loadModel(par, model);
	if (algo == 1)
		par->relaxIntegrality_[1] = true;

	/** solve */
	solver->solve();

#ifdef USE_MPI
	if (algo != 2 || MPI::COMM_WORLD.Get_rank() == 0)
#endif
	{
		cout << "Total solution time: " << solver->solutionTime_ << " seconds" << endl;
		cout << "Solution status: " << solver->status_ << endl;
		if (algo != 2)
			cout << "Objective value: " << solver->primalBound_ << endl;
		else
		{
			if (feasibleRecovery >= 0) cout << "Best upper bound: " << solver->primalBound_ << endl;
			cout << "Best lower bound: " << solver->dualBound_ << endl;
		}
	}

	/** free memory */
	if (solver)
	{
		delete solver;
		solver = NULL;
	}
	if (par)
	{
		delete par;
		par = NULL;
	}
	if (model)
	{
		delete model;
		model = NULL;
	}

#ifdef USE_MPI
	if (algo == 2) MPI_Finalize();
#endif

	return 0;
}
