/*
 * dsp.cpp
 *
 *  Created on: Feb 21, 2018
 *      Author: Kibaek Kim
 */

#include <iostream>
#include "CoinMpsIO.hpp"
#include "DspCInterface.h"

const char* gDspUsage = 
	"Not enough or invalid arguments, please try again.\n\n"
	"Usage: --algo <de,bd,dd,drbd,drdd,dw> --smps <smps file> --mps <mps file> --dec <dec file> [--soln <solution file prefix> --param <param file> --test <benchmark objective value>]\n\n"
	"       --algo\tchoice of algorithms.\n"
	"             \tde: deterministic equivalent form\n"
	"             \tbd: Benders decomposition\n"
	"             \tdd: dual decomposition\n"
	"             \tdrbd: distributionally robust bd\n"
	"             \tdrdd: distributionally robust dd\n"
	"             \tdw: Dantzig-Wolfe decomposition with branch-and-bound\n"
	"       --smps\tSMPS file name without extensions. For example, if your SMPS files are ../test/farmer.cor, ../test/farmer.sto, and ../test/farmer.tim, this value should be ../test/farmer\n"
	"       --mps\tMPS file name\n"
	"       --dec\tDEC file name\n"
	"       --soln\toptional argument for solution file prefix. For example, if the prefix is given as MySol, then two files MySol.primal.txt and MySol.dual.txt will be written for primal and dual solutions, respectively.\n"
	"       --param\toptional paramater for parameter file name\n"
	"       --test\toptional parameter for testing objective value\n";

void setBlockIds(DspApiEnv* env, int nsubprobs, bool master_has_subblocks);
int runDsp(char* algotype, char* smpsfile, char* mpsfile, char* decfile, char* solnfile, char* paramfile, char* testvalue);
int readMpsDec(DspApiEnv* env, char* mpsfile, char* decfile);
int parseDecFile(char* decfile, vector<vector<string> >& rows_in_blocks);
void createBlockModel(DspApiEnv* env, CoinMpsIO& p, const CoinPackedMatrix* mat, 
	int blockid, vector<string>& rows_in_block, map<string,int>& rowname2index, 
	const char* ctype, const double* obj);

const double test_tolerance = 1.0e-2;

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
		char* mpsfile = NULL;
		char* decfile = NULL;
		char* solnfile = NULL;
		char* paramfile = NULL;
		char* testvalue = NULL;
		for (int i = 1; i < argc; i += 2) {
			if (i + 1 != argc) {
				if (string(argv[i]) == "--algo") {
					algotype = argv[i+1];
				} else if (string(argv[i]) == "--smps") {
					smpsfile = argv[i+1];
				} else if (string(argv[i]) == "--mps") {
					mpsfile = argv[i+1];
				} else if (string(argv[i]) == "--dec") {
					decfile = argv[i+1];
				} else if (string(argv[i]) == "--soln") {
					solnfile = argv[i+1];
				} else if (string(argv[i]) == "--param") {
					paramfile = argv[i+1];
				} else if (string(argv[i]) == "--test") {
					testvalue = argv[i+1];
				} else {
					EXIT_WITH_MSG
				}
			}
		}

		// algotype is required.
		if (algotype == NULL) {
			EXIT_WITH_MSG
		}

		// Either smps or mps/dec files are required.
		if (smpsfile == NULL && (mpsfile == NULL || decfile == NULL)) {
			EXIT_WITH_MSG
		}

		// run dsp
		int ret = runDsp(algotype, smpsfile, mpsfile, decfile, solnfile, paramfile, testvalue);

#ifdef DSP_HAS_MPI
		MPI_Finalize();
#endif
		return ret;
	}
#undef EXIT_WITH_MSG
}

int runDsp(char* algotype, char* smpsfile, char* mpsfile, char* decfile, char* solnfile, char* paramfile, char* testvalue) {

	int ret = 0;
	bool isroot = true;
	bool issolved = true;
	bool isstochastic = true;
#ifdef DSP_HAS_MPI
	int comm_rank, comm_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
	isroot = comm_rank == 0 ? true : false;
#endif

	if (isroot) cout << "Creating DSP environment\n";
	DspApiEnv* env = createEnv();

	// Read problem instance from file(s)
	if (smpsfile != NULL) {
		if (isroot) cout << "Reading SMPS files: " << smpsfile << endl;
		ret = readSmps(env, smpsfile);
		if (ret != 0) return ret;
		if (isroot) {
			cout << "First stage: " << getNumRows(env,0) << " rows, " << getNumCols(env,0) << " cols, " << getNumIntegers(env,0) << " integers" << endl;
			cout << "Second stage: " << getNumRows(env,1) << " rows, " << getNumCols(env,1) << " cols, " << getNumIntegers(env,1) << " integers" << endl;
			cout << "Number of scenarios: " << getNumSubproblems(env) << endl;
		}
		setBlockIds(env, getNumSubproblems(env), true);
	} else if (mpsfile != NULL && decfile != NULL) {
		if (isroot) {
			cout << "Reading MPS file: " << mpsfile << endl;
			cout << "Reading DEC file: " << decfile << endl;
		}
		ret = readMpsDec(env, mpsfile, decfile);
		if (ret != 0) return ret;
		isstochastic = false;
	}

	if (paramfile != NULL) {
		if (isroot) cout << "Reading parameter files: " << paramfile << endl;
		readParamFile(env, paramfile);
	}

	if (string(algotype) == "de") {
		solveDe(env);
	} else if (string(algotype) == "bd") {
		if (isstochastic) {
#ifdef DSP_HAS_MPI
			solveBdMpi(env, MPI_COMM_WORLD);
#else
			solveBd(env);
#endif
		} else {
			cout << "Benders decomposition is not available for mps/dec files." << endl;
			issolved = false;
		}
	} else if (string(algotype) == "dd") {
		if (isstochastic) {
#ifdef DSP_HAS_MPI
			solveDdMpi(env, MPI_COMM_WORLD);
#else
			solveDd(env);
#endif
		} else {
			cout << "Dual decomposition is not available for mps/dec files." << endl;
			issolved = false;
		}
	}
	else if (string(algotype) == "drbd")
	{
		if (isstochastic) {
			char drofile[128];
			sprintf(drofile, "%s.dro", smpsfile);
			if (isroot) cout << "Reading DRO file: " << drofile << endl;

			int ret = readDro(env, drofile);
			if (ret != 0) return ret;

			if (isroot)
				cout << "Run distributionally robust Benders decomposition" << endl;
#ifdef DSP_HAS_MPI
			solveBdMpi(env, MPI_COMM_WORLD);
#else
			solveBd(env);
#endif
		} else {
			cout << "Benders decomposition is not available for mps/dec files." << endl;
			issolved = false;
		}
	}
	else if (string(algotype) == "drdd")
	{
		if (isstochastic) {
			char drofile[128];
			sprintf(drofile, "%s.dro", smpsfile);
			if (isroot) cout << "Reading DRO file: " << drofile << endl;

			int ret = readDro(env, drofile);
			if (ret != 0) return ret;

			if (isroot)
				cout << "Run distributionally robust dual decomposition" << endl;
			solveDd(env);
		} else {
			cout << "Dual decomposition is not available for mps/dec files." << endl;
			issolved = false;
		}
	}
	else if (string(algotype) == "dw")
	{
#ifdef DSP_HAS_MPI
		solveDwMpi(env, MPI_COMM_WORLD);
#else
		solveDw(env);
#endif
	}
	else
	{
		if (isroot) cout << "Invalid algorithm type, please try again.\n";
		issolved = false;
	}

	/** parse result */
	if (isroot && issolved) {
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
			cout << "Gap (%)     : " << fabs(primobj-dualobj)/(fabs(primobj)+1.e-10)*100 << endl;

			if (testvalue != NULL) {
				double val = atof(testvalue);
				cout << "Testing Bound: " << val << endl;
				// cout << "relgap: " << fabs(val - primobj) / (fabs(val) + 1.e-10) << endl;
				// cout << "relgap: " << fabs(val - dualobj) / (fabs(val) + 1.e-10) << endl;
				if (primobj >= dualobj) {
					if ((val - primobj) / (fabs(val) + 1.e-10) > test_tolerance || (dualobj - val) / (fabs(val) + 1.e-10) > test_tolerance)
						ret = 1;
				} else {
					if ((primobj - val) / (fabs(val) + 1.e-10) > test_tolerance || (val - dualobj) / (fabs(val) + 1.e-10) + test_tolerance)
						ret = 1;
				}
			}

			/** write solutions to files */
			if (solnfile != NULL) {

				char pobjname[128];
				sprintf(pobjname, "%s.primobj.txt", solnfile);
				ofstream pobjstream(pobjname);
				pobjstream << primobj << endl;
				pobjstream.close();

				char dobjname[128];
				sprintf(dobjname, "%s.dualobj.txt", solnfile);
				ofstream dobjstream(dobjname);
				dobjstream << dualobj << endl;
				dobjstream.close();

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

	return ret;
}

void setBlockIds(DspApiEnv* env, int nsubprobs, bool master_has_subblocks) {
	int comm_rank = 0, comm_size = 1;
#ifdef DSP_HAS_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
#endif

	// empty block IDs
	vector<int> proc_idx_set;

	// DSP is further parallelized with comm_size > nsubprobs.
	int modrank = comm_rank % nsubprobs;
	if (master_has_subblocks == false) {
		if (comm_size == 1) {
			// assign sub-blocks in round-robin fashion
			for (int s = modrank; s < nsubprobs; s += comm_size)
				proc_idx_set.push_back(s);
		} else {
			if (comm_rank > 0) {
				// exclude master
				comm_size--;
				modrank = (comm_rank-1) % nsubprobs;
	
				// assign sub-blocks in round-robin fashion
				for (int s = modrank; s < nsubprobs; s += comm_size) {
					proc_idx_set.push_back(s);
				}
			}
		}
	} else {
		// assign sub-blocks in round-robin fashion
		for (int s = modrank; s < nsubprobs; s += comm_size) {
			proc_idx_set.push_back(s);
		}
#if 0
		for (int i = 0; i < comm_size; ++i) {
			if (i == comm_rank) {
				printf("comm_rank %d:\n", i);
				for (int j = 0; j < proc_idx_set.size(); ++j)
					printf("  %d\n", proc_idx_set[j]);
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
#endif
	}

	// set the block ids to Dsp
	setIntPtrParam(env, "ARR_PROC_IDX", (int) proc_idx_set.size(), &proc_idx_set[0]);
}

int readMpsDec(DspApiEnv* env, char* mpsfile, char* decfile) {
	int ret = 0;
	// Read .mps file
	CoinMpsIO p;
	p.readMps(mpsfile);
	const CoinPackedMatrix* mps_matrix = p.getMatrixByRow();
	int ncols = p.getNumCols();
	const double* obj = p.getObjCoefficients();
	vector<double> zeros(ncols, 0.0);
	string ctype;
	for (int j = 0; j < ncols; ++j)
		if (p.isContinuous(j))
			ctype.push_back('C');
		else
			ctype.push_back('I');

	//cout << "Creating hash table for matrix rows ... ";
	map<string, int> rowname2index;
	for (int i = 0; i < p.getNumRows(); ++i)
		rowname2index[string(p.rowName(i))] = i;
	//cout << "done!" << endl;

	// Read .dec file
	// For each block, the following vector stores the corresponding rows.
	//cout << "Parsing .dec file ... ";
	vector<vector<string> > rows_in_blocks;
	ret = parseDecFile(decfile, rows_in_blocks);
	if (ret != 0) return ret;

	// Assign block(s) to each process
	setBlockIds(env, rows_in_blocks.size() - 1, true);

	vector<int> proc_idx_set;
	for (int s = 0; s < env->par_->getIntPtrParamSize("ARR_PROC_IDX"); ++s) {
		proc_idx_set.push_back(env->par_->getIntPtrParam("ARR_PROC_IDX")[s]);
	}

	// For master block
	int blockid = 0;
	createBlockModel(env, p, mps_matrix, 0, rows_in_blocks[0], rowname2index, ctype.c_str(), obj);

	// For each sub-block
	for (unsigned i = 0; i < proc_idx_set.size(); ++i) {
		blockid = proc_idx_set[i] + 1; // blockid = subproblem index + 1, because master = 0
		createBlockModel(env, p, mps_matrix, blockid, rows_in_blocks[blockid], rowname2index, ctype.c_str(), obj);
	}

	updateBlocks(env);

	return ret;
}

int parseDecFile(char* decfile, vector<vector<string> >& rows_in_blocks) {
	// Cards defined in .dec file
	enum DecCard {
		PRESOLVED = 0,
		NBLOCKS,
		BLOCK,
		MASTERCONSS
	} card;

	int nblocks, current_block = -1;
	rows_in_blocks.clear();
	string line;
	ifstream myfile(decfile);
	if (myfile.is_open()) {
		string card_prefix;
		size_t found;
		while(getline(myfile, line)) {
			found = line.find_first_of(" ");
			card_prefix = line.substr(0, found);
			if (card_prefix.compare("PRESOLVED") == 0)
				card = PRESOLVED;
			else if (card_prefix.compare("NBLOCKS") == 0)
				card = NBLOCKS;
			else if (card_prefix.compare("BLOCK") == 0) {
				card = BLOCK;
				current_block++;
			} else if (card_prefix.compare("MASTERCONSS") == 0)
				card = MASTERCONSS;
			else {
				if (card == NBLOCKS) {
					nblocks = atoi(line.c_str());
					rows_in_blocks.resize(nblocks+1);
				} else if (card == BLOCK)
					rows_in_blocks[current_block+1].push_back(line);
				else if (card == MASTERCONSS)
					rows_in_blocks[0].push_back(line);
			}
		}
	} else {
		printf("Cannot open file: %s\n", decfile);
		return 1;
	}
#if 0
	cout << "Number of blocks: " << rows_in_blocks.size() << endl;
	for (int i = 0; i < rows_in_blocks.size(); ++i) {
		cout << "Block " << i << endl;
		for (int j = 0; j < rows_in_blocks[i].size(); ++j)
			cout << rows_in_blocks[i][j] << endl;
	}
#endif
	return 0;
}

void createBlockModel(DspApiEnv* env, CoinMpsIO& p, const CoinPackedMatrix* mat, 
		int blockid, vector<string>& rows_in_block, map<string,int>& rowname2index, 
		const char* ctype, const double* obj) {
	vector<int> rowids(rows_in_block.size(),-1);
	vector<double> rlbd(rows_in_block.size(),0.0);
	vector<double> rubd(rows_in_block.size(),0.0);
	CoinPackedMatrix submat(false,0,0);

	//cout << "Creating block " << blockid << " ... ";
	for (unsigned j = 0; j < rows_in_block.size(); ++j) {
		int k = rowname2index[rows_in_block[j]];
		rowids[j] = k;
		rlbd[j] = p.getRowLower()[k];
		rubd[j] = p.getRowUpper()[k];
	}
	//cout << " with " << rowids.size() << " rows ... ";
	submat.submatrixOf(*mat, rowids.size(), &rowids[0]);
	//cout << "done!" << endl;
#if 0
	CoinMpsIO pout;
	vector<string> cnames;
	for (int j = 0; j < p.getNumCols(); ++j)
		cnames.push_back(string(p.columnName(j)));
	pout.setMpsData(submat, 1.0e+20, p.getColLower(), p.getColUpper(), obj, ctype, &rlbd[0], &rubd[0], cnames, rows_in_block);
	char pout_name[128];
	sprintf(pout_name, "block%d.mps", blockid);
	pout.writeMps(pout_name);
#endif

	//cout << "Loading the block ... ";
	loadBlockProblem(env, blockid, p.getNumCols(), rowids.size(), 
		submat.getNumElements(), submat.getVectorStarts(), submat.getIndices(), submat.getElements(), 
		p.getColLower(), p.getColUpper(), ctype, obj, &rlbd[0], &rubd[0]);
	//cout << "done!" << endl;
}

