/*
 * DspParams.cpp
 *
 *  Created on: Oct 20, 2015
 *      Author: kibaekkim
 */

//#define DSP_DEBUG
#include "Utility/DspMessage.h"
#include "Utility/DspParams.h"

#define MAX_INT_NUM numeric_limits<int>::max()
#define MAX_DBL_NUM numeric_limits<double>::max()

DspParams::DspParams()
{
	initBoolParams();
	initIntParams();
	initDblParams();
	initStrParams();
	initBoolPtrParams();
	initIntPtrParams();
}

DspParams::~DspParams()
{
	BoolPtrParams_.deleteParam("RELAX_INTEGRALITY");
	IntPtrParams_.deleteParam("BD/ARR_AUG_SCENS");
	IntPtrParams_.deleteParam("ARR_PROC_IDX");
}

/** read parameter file */
void DspParams::readParamFile(const char * param_file)
{
	string line;
	ifstream myfile(param_file);
	if (myfile.is_open())
	{
		string param_element[3];
		size_t startpos, found;
		while(getline(myfile, line))
		{
			//printf("Read line: %s\n", line.c_str());
			/** comment out? */
			startpos = line.find_first_of("#");
			if(string::npos != startpos) line = line.substr(0, startpos);

			/** empty line? */
			startpos = line.find_first_not_of(" \t");
			if(string::npos == startpos) continue;

			bool is_valid = true;
			for (int i = 0; i < 3; ++i)
			{
				/** ltrim */
				startpos = line.find_first_not_of(" \t");
				if(string::npos != startpos) line = line.substr(startpos);

				/** read element */
				found = line.find_first_of(" \t");
				if(string::npos == found && i < 2)
				{
					printf("Invalid parameter format.\n");
					is_valid = false;
					break;
				}
				param_element[i] = line.substr(0, found);
				if (string::npos != found)
					line = line.substr(found);
			}
			if (!is_valid) break;

			if (param_element[0].compare("bool") == 0)
			{
				if (param_element[2].compare("true") == 0)
					BoolParams_.setParam(param_element[1], true);
				else if (param_element[2].compare("false") == 0)
					BoolParams_.setParam(param_element[1], false);
				else
				{
					printf("Invalid parameter type.\n");
					is_valid = false;
					break;
				}
			}
			else if (param_element[0].compare("double") == 0)
			{
				DblParams_.setParam(param_element[1], atof(param_element[2].c_str()));
			}
			else if (param_element[0].compare("int") == 0)
			{
				IntParams_.setParam(param_element[1], atoi(param_element[2].c_str()));
			}
			else if (param_element[0].compare("string") == 0)
			{
				StrParams_.setParam(param_element[1], param_element[2]);
			}
			else
			{
				printf("Invalid parameter type.\n");
				is_valid = false;
				break;
			}
		}
		myfile.close();
	}
	else
		printf("Unable to open parameter file <%s>.\n", param_file);
}

void DspParams::initBoolParams()
{
	/** enable trust region */
	BoolParams_.createParam("DD/TR", true);

	/** enable decreasing trust region */
	BoolParams_.createParam("DD/TR/DECREASE", true);

	/** enable decreasing trust region */
	BoolParams_.createParam("DD/ALLOW_IDLE_WORKERS", true);

	/** cache recourse models */
	BoolParams_.createParam("DD/CACHE_RECOURSE", true);

	/** log dual variable values */
	BoolParams_.createParam("DD/LOG_DUAL_VARS", false);

	/** log dual variable values */
	BoolParams_.createParam("DD/LOG_LB_TIME", false);

	/** enable asynchronous parallelization */
	BoolParams_.createParam("DD/ASYNC", false);

	/** dynamic allocation in the asynchronous DD */
	BoolParams_.createParam("DD/ASYNC/DYNAMIC", false);

	/** static FIFO scheduling in the asynchronous DD; otherwise LIFO */
	BoolParams_.createParam("DD/ASYNC/FIFO", true);

	/** options for Dantzig-Wolfe decomposition */
	BoolParams_.createParam("DW/MASTER/PIPS", false);
	BoolParams_.createParam("DW/MASTER/IPM", false);
	BoolParams_.createParam("DW/MASTER/BRANCH_ROWS", false);
	BoolParams_.createParam("DW/MASTER/REUSE_COLS", false);
	BoolParams_.createParam("DW/TRUST_REGION", false);
	BoolParams_.createParam("DW/HEURISTICS", true);
	BoolParams_.createParam("DW/HEURISTICS/ROUNDING", false);
	BoolParams_.createParam("DW/HEURISTICS/SMIP", true);
	BoolParams_.createParam("DW/STRONG_BRANCH", false);
	BoolParams_.createParam("DW/BRANCH/INTEGER_FIRST", false);
}

void DspParams::initIntParams()
{
	/** print level */
	IntParams_.createParam("LOG_LEVEL", 1);

	/** branch-and-cut node limit */
	IntParams_.createParam("NODE_LIM", MAX_INT_NUM);

	/** number of cuts to the master per iteration */
	IntParams_.createParam("BD/NUM_CUTS_PER_ITER", 1);

	/** number of cores used in OpenMP library (Benders only) */
	IntParams_.createParam("NUM_CORES", 1);

	/** Benders cut priority (refer CONSHDLR_SEPAPRIORITY of SCIP constraint handler */
	IntParams_.createParam("BD/CUT_PRIORITY", -200000);

	/** Benders lower bound methods:
	 * 0 = solve separate LP relaxation problems;
	 * 1 = solve separate MILP relaxation problems */
	IntParams_.createParam("BD/INIT_LB_ALGO", SEPARATE_LP);

    /** iteration limit of the dual decomposition used in Benders for initial lower bounding */
	IntParams_.createParam("BD/DD/ITER_LIM", 1);

	/** iteration limit */
    IntParams_.createParam("DD/ITER_LIM", MAX_INT_NUM);

	/** algorithm for the master */
#ifdef DSP_HAS_OOQP
	IntParams_.createParam("DD/MASTER_ALGO", IPM_Feasible);
#else
	IntParams_.createParam("DD/MASTER_ALGO", IPM);
#endif

	/** number of cuts to the master per iteration */
	IntParams_.createParam("DD/NUM_CUTS_PER_ITER", 1);

	/** add feasibility cuts */
	IntParams_.createParam("DD/FEAS_CUTS", 1);

	/** add optimality cuts */
	IntParams_.createParam("DD/OPT_CUTS", 1);

	/** evaluate upper bound */
	IntParams_.createParam("DD/EVAL_UB", 1);
	IntParams_.createParam("DW/EVAL_UB", 0);

	/** maximum number of solutions to evaluate */
	IntParams_.createParam("DD/MAX_EVAL_UB", 1);
	IntParams_.createParam("DW/MAX_EVAL_UB", 1);

	/** maximum queue size for asynchronous one */
	IntParams_.createParam("DD/MAX_QSIZE", 5);

	/** minimum number of processes to wait at the master */
	IntParams_.createParam("DD/MIN_PROCS", 1);

//#ifdef DSP_HAS_GRB
//	IntParams_.createParam("DE/SOLVER", OsiGrb);
	
//#endif

#ifdef DSP_HAS_CPX
	IntParams_.createParam("BD/SUB/SOLVER", OsiCpx);
	IntParams_.createParam("DD/MASTER/SOLVER", OsiCpx);
	IntParams_.createParam("DD/SUB/SOLVER", OsiCpx);
	IntParams_.createParam("DE/SOLVER", OsiCpx);
	IntParams_.createParam("DW/MASTER/SOLVER", OsiCpx);
	IntParams_.createParam("DW/SUB/SOLVER", OsiCpx);
#else
#ifdef DSP_HAS_GRB
	IntParams_.createParam("BD/SUB/SOLVER", OsiGrb);
	IntParams_.createParam("DE/SOLVER", OsiGrb);
	IntParams_.createParam("DD/MASTER/SOLVER", OsiGrb);
	IntParams_.createParam("DD/SUB/SOLVER", OsiGrb);
	IntParams_.createParam("DW/MASTER/SOLVER", OsiGrb);
	IntParams_.createParam("DW/SUB/SOLVER", OsiGrb);
#else
	IntParams_.createParam("DW/MASTER/SOLVER", OsiClp);
	IntParams_.createParam("BD/SUB/SOLVER", OsiClp);
	IntParams_.createParam("DD/MASTER/SOLVER", OsiClp);
#ifdef DSP_HAS_SCIP
	IntParams_.createParam("DD/SUB/SOLVER", OsiScip);
	IntParams_.createParam("DE/SOLVER", OsiScip);
	IntParams_.createParam("DW/SUB/SOLVER", OsiScip);
#endif
#endif
#endif


	IntParams_.createParam("DE/SOLVER/LOG_LEVEL", 0);
	IntParams_.createParam("DD/SUB/SOLVER/LOG_LEVEL", 0);
	IntParams_.createParam("DW/MASTER/SOLVER/LOG_LEVEL", 0);
	IntParams_.createParam("DW/SUB/SOLVER/LOG_LEVEL", 0);

	/** number of threads used for subproblem solution */
	IntParams_.createParam("BD/MASTER/THREADS", 1);
	IntParams_.createParam("BD/SUB/THREADS", 1);
	IntParams_.createParam("DD/SUB/THREADS", 1);
	IntParams_.createParam("DW/SUB/THREADS", 1);

	/** display frequency */
	IntParams_.createParam("SCIP/DISPLAY_FREQ", 100);

	IntParams_.createParam("ALPS/SEARCH_STRATEGY", 0);
	IntParams_.createParam("ALPS/NODE_LIM", 1000000);
	IntParams_.createParam("ALPS/NODE_LOG_INTERVAL", 1);

	IntParams_.createParam("DW/MASTER/COL_AGE_LIM", 10);
    IntParams_.createParam("DW/ITER_LIM", MAX_INT_NUM);
    IntParams_.createParam("DW/HEURISTICS/TRIVIAL/ITER_LIM", MAX_INT_NUM);
	IntParams_.createParam("DW/HEURISTICS/DIVE/ITER_LIM", MAX_INT_NUM);
	IntParams_.createParam("DW/SUB/ADVIND", 1);
	IntParams_.createParam("DW/BRANCH", 2);
	IntParams_.createParam("DW/STRONG_BRANCH/ITER_LIM", 10);
}

void DspParams::initDblParams()
{
	/** wall clock limit */
	DblParams_.createParam("DE/WALL_LIM", MAX_DBL_NUM);

	/** wall clock limit */
	DblParams_.createParam("BD/WALL_LIM", MAX_DBL_NUM);

	/** wall clock limit */
	DblParams_.createParam("DD/WALL_LIM", MAX_DBL_NUM);

	/** TODO: is this wall clock limit? */
	DblParams_.createParam("DW/TIME_LIM", MAX_DBL_NUM);

	/** options for trust region */
	DblParams_.createParam("DD/TR/SIZE", 100.0);
	DblParams_.createParam("DW/TR/SIZE", 100.0);
	DblParams_.createParam("DW/MIN_INCREASE", 1.0e-6);
	DblParams_.createParam("DW/INIT_CENTER", 0.1);

	/** stopping tolerance */
	DblParams_.createParam("DD/STOP_TOL", 0.00001);

	/** branch-and-bound gap tolerance */
	DblParams_.createParam("DE/GAPTOL", 1.0e-4);
	DblParams_.createParam("DD/SUB/GAPTOL", 1.0e-5);
	DblParams_.createParam("DW/GAPTOL", 1.0e-4);
	DblParams_.createParam("DW/SUB/GAPTOL", 1.0e-4);

	/** time limit */
	DblParams_.createParam("DD/SUB/TIME_LIM", 1e+20);
	DblParams_.createParam("DW/SUB/TIME_LIM", 1e+20);

	/** LB-UB worker ratio */
	DblParams_.createParam("DD/WORKER_RATIO", 0.8);

	/** minimum wait time for the master to receive worker processes */
	DblParams_.createParam("DD/ASYNC/MIN_WAIT_TIME", 5.0);

	/** options for branch-and-bound search */
	DblParams_.createParam("ALPS/TIME_LIM", MAX_DBL_NUM);
	DblParams_.createParam("DW/HEURISTICS/TRIVIAL/TIME_LIM", MAX_DBL_NUM);
	DblParams_.createParam("DW/HEURISTICS/DIVE/TIME_LIM", MAX_DBL_NUM);
}

void DspParams::initStrParams()
{
	/** prefix for output files */
	StrParams_.createParam("OUTPUT/PREFIX", "dsp");
	StrParams_.createParam("DW/LOGFILE/OBJS", "");
	StrParams_.createParam("VBC/FILE", "");
}

void DspParams::initBoolPtrParams()
{
	/** array of indicating relaxation of integer variables;
	 * each element represents a stage. */
	BoolPtrParams_.createParam("RELAX_INTEGRALITY", 2);
	BoolPtrParams_.setParam("RELAX_INTEGRALITY", 0, false);
	BoolPtrParams_.setParam("RELAX_INTEGRALITY", 1, false);
}

void DspParams::initIntPtrParams()
{
	DSPdebugMessage("creating integer array parameters\n");
	/** array of scenarios for the current process */
	IntPtrParams_.createParam("ARR_PROC_IDX");

	/** array of augmented scenarios */
	IntPtrParams_.createParam("BD/ARR_AUG_SCENS");
}
