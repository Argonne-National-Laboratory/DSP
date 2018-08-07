/*
 * dsp.cpp
 *
 *  Created on: Feb 21, 2018
 *      Author: Kibaek Kim
 */

#include <iostream>
#ifdef DSP_HAS_MPI
#include <mpi.h>
#endif
#include <DspCInterface.h>

const char* gDspUsage = 
	"Not enough or invalid arguments, please try again.\n\n"
	"Usage: --algo <de,bd,dd> --smps <smps file> [--soln <solution file prefix> --param <param file>]\n\n"
	"       --algo\tchoice of algorithms. de: deterministic equivalent form; bd: Benders decomposition; dd: dual decomposition\n"
	"       --smps\tSMPS file name without extensions. For example, if your SMPS files are ../test/farmer.cor, ../test/farmer.sto, and ../test/farmer.tim, this value should be ../test/farmer\n"
	"       --soln\toptional argument for solution file prefix. For exampe, if the prefix is given as MySol, then two files MySol.primal.txt and MySol.dual.txt will be written for primal and dual solutions, respectively.\n"
	"       --param\toptional paramater for parameter file name\n";

void setBlockIds(DspApiEnv* env, bool master_has_subblocks = false);
void runDsp(char* algotype, char* smpsfile, char* solnfile, char* paramfile);

/*
 This will compile a stand-alone binary file that reads problem instances.
*/
int main(int argc, char* argv[]) {

	bool isroot = true;
#ifdef DSP_HAS_MPI
	int comm_rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
	isroot = comm_rank == 0 ? true : false;
#define EXIT_WITH_MSG \
	if(isroot) cout << gDspUsage; \
	MPI_Finalize(); \
	exit(0); 
#else
#define EXIT_WITH_MSG \
	if(isroot) cout << gDspUsage; \
	exit(0); 
#endif


	if (argc < 5) {
		EXIT_WITH_MSG
	} else {
		char* algotype = NULL;
		char* smpsfile = NULL;
		char* solnfile = NULL;
		char* paramfile = NULL;
		for (int i = 1; i < argc; ++i) {
			if (i + 1 != argc) {
				if (string(argv[i]) == "--algo") {
					algotype = argv[i+1];
				} else if (string(argv[i]) == "--smps") {
					smpsfile = argv[i+1];
				} else if (string(argv[i]) == "--soln") {
					solnfile = argv[i+1];
				} else if (string(argv[i]) == "--param") {
					paramfile = argv[i+1];
				} else {
					EXIT_WITH_MSG
				}
			}
			i++;
		}

		if (algotype == NULL || smpsfile == NULL) {
			EXIT_WITH_MSG
		}

		/** run dsp */
		runDsp(algotype, smpsfile, solnfile, paramfile);

#ifdef DSP_HAS_MPI
		MPI_Finalize();
#endif

		return 0;
	}
#undef EXIT_WITH_MSG
}

void runDsp(char* algotype, char* smpsfile, char* solnfile, char* paramfile) {

	bool isroot = true;
#ifdef DSP_HAS_MPI
	int comm_rank, comm_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
	isroot = comm_rank == 0 ? true : false;
#endif

	if (isroot) cout << "Creating DSP environment\n";
	DspApiEnv* env = createEnv();

	if (isroot) cout << "Reading SMPS files: " << smpsfile << endl;
	int ret = readSmps(env, smpsfile);
	if (ret != 0) return;
	setBlockIds(env, false);

	if (paramfile != NULL) {
		if (isroot) cout << "Reading parameter files: " << paramfile << endl;
		readParamFile(env, paramfile);
	}

	if (string(algotype) == "de")
		solveDe(env);
#ifdef DSP_HAS_MPI
	else if (string(algotype) == "bd") {
		if (comm_size > 1)
			solveBdMpi(env, MPI_COMM_WORLD);
		else
			solveBd(env);
	} else if (string(algotype) == "dd") {
		if (comm_size > 1)
			solveDdMpi(env, MPI_COMM_WORLD);
		else
			solveDd(env);
	}
#else
	else if (string(algotype) == "bd")
		solveBd(env);
	else if (string(algotype) == "dd")
		solveDd(env);
#endif
	else {
		if (isroot) cout << "Invalid algorithm type, please try again.\n";
	}

	/** parse result */
	if (isroot) {
		int status = getStatus(env);
		cout << "Status: " << status << endl;
		if (status == DSP_STAT_OPTIMAL ||
			status == DSP_STAT_LIM_ITERorTIME ||
			status == DSP_STAT_STOPPED_GAP) {

			double primobj = getPrimalBound(env);
			double dualobj = getDualBound(env);
			int nrows = getNumCouplingRows(env);
			int ncols = getTotalNumCols(env);

			cout << "Primal Bound: " << primobj << endl;
			cout << "Dual Bound  : " << dualobj << endl;

			/** write solutions to files */
			if (solnfile != NULL) {
				if (primobj < 1.0e+20) {
					double *primsol = new double [ncols];
					getPrimalSolution(env, ncols, primsol);

					char psolname[128];
					sprintf(psolname, "%s.primal.txt", solnfile);
					ofstream psolstream(psolname);
					//cout << "Priaml Solution:" << endl;
					for (int j = 0; j < ncols; ++j)
						psolstream << primsol[j] << endl;
					psolstream.close();

					delete [] primsol;
				}
	
				if (string(algotype) == "dd") {
					double *dualsol = new double [nrows];
					getDualSolution(env, nrows, dualsol);
	
					char dsolname[128];
					sprintf(dsolname, "%s.dual.txt", solnfile);
					ofstream dsolstream(dsolname);
					//cout << "Dual Solution:" << endl;
					for (int j = 0; j < nrows; ++j)
						dsolstream << dualsol[j] << endl;
					dsolstream.close();
					
					delete [] dualsol;
				}
			}
		}
	}


	if (isroot) cout << "Deleting DSP environment\n";
	freeEnv(env);
}

void setBlockIds(DspApiEnv* env, bool master_has_subblocks) {
#ifdef DSP_HAS_MPI
	int nsubprobs = getNumSubproblems(env);
	int comm_rank, comm_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

	// empty block IDs
	vector<int> proc_idx_set;

	// DSP is further parallelized with comm_size > nsubprobs.
	int modrank = comm_rank % nsubprobs;
	// If we have more than one processor, do not assign a sub-block to the master.
	if (master_has_subblocks == false && comm_rank > 0) {
		// exclude master
		comm_size--;
		modrank = (comm_rank-1) % nsubprobs;
	}

	// assign sub-blocks in round-robin fashion
	for (int s = modrank; s < nsubprobs; s += comm_size)
		proc_idx_set.push_back(s);

	// set the block ids to Dsp
	setIntPtrParam(env, "ARR_PROC_IDX", (int) proc_idx_set.size(), &proc_idx_set[0]);
#endif
}

