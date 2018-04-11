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
#include <CoinMpsIO.hpp>
#include <DspCInterface.h>

const char* gDspUsage = 
	"Not enough or invalid arguments, please try again.\n\n"
	"Usage: --smps <smps file> --mps <mps file> --dec <dec file> [--param <param file>]\n\n"
	"       --smps\tSMPS file name without extensions. For example, if your SMPS files are ../test/farmer.cor, ../test/farmer.sto, and ../test/farmer.tim, this value should be ../test/farmer\n"
	"       --mps\tMPS file name\n"
	"       --dec\tDEC file name\n"
	"       --param\toptional paramater for parameter file name\n";

void setBlockIds(DspApiEnv* env, int nblocks, bool master_has_subblocks = false);
void runDsp(char* smpsfile, char* mpsfile, char* decfile, char* paramfile);
void readMpsDec(DspApiEnv* env, char* mpsfile, char* decfile);
void parseDecFile(char* decfile, vector<vector<string> >& rows_in_blocks);
void createBlockModel(DspApiEnv* env, CoinMpsIO& p, const CoinPackedMatrix* mat, 
	int blockid, vector<string>& rows_in_block, const char* ctype, const double* obj);

std::vector<int> g_proc_idx_set;

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
	if(isroot) std::cout << gDspUsage; \
	MPI_Finalize(); \
	exit(0); 
#else
#define EXIT_WITH_MSG \
	if(isroot) std::cout << gDspUsage; \
	exit(0); 
#endif

	if (argc < 5) {
		EXIT_WITH_MSG
	} else {
		char* smpsfile = NULL;
		char* mpsfile = NULL;
		char* decfile = NULL;
		char* paramfile = NULL;
		for (int i = 1; i < argc; ++i) {
			if (i + 1 != argc) {
				if (std::string(argv[i]) == "--smps") {
					smpsfile = argv[i+1];
				} else if (std::string(argv[i]) == "--mps") {
					mpsfile = argv[i+1];
				} else if (std::string(argv[i]) == "--dec") {
					decfile = argv[i+1];
				} else if (std::string(argv[i]) == "--param") {
					paramfile = argv[i+1];
				} else {
					EXIT_WITH_MSG
				}
			}
			i++;
		}

		if (smpsfile == NULL && (mpsfile == NULL || decfile == NULL)) {
			EXIT_WITH_MSG
		}

		/** run dsp */
		runDsp(smpsfile, mpsfile, decfile, paramfile);

#ifdef DSP_HAS_MPI
		MPI_Finalize();
#endif

		return 0;
	}
#undef EXIT_WITH_MSG
}

void runDsp(char* smpsfile, char* mpsfile, char* decfile, char* paramfile) {
	bool isroot = true;
#ifdef DSP_HAS_MPI
	int comm_rank, comm_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
	isroot = comm_rank == 0 ? true : false;
#endif

	if (isroot) cout << "Creating DSP environment" << endl;
#ifdef DSP_HAS_MPI
	DspApiEnv* env = createEnv();
#endif

	// Read problem instance from file(s)
	if (smpsfile != NULL) {
		if (isroot) cout << "Reading SMPS files: " << smpsfile << endl;
		readSmps(env, smpsfile);
		setBlockIds(env, getNumScenarios(env));
	} else if (mpsfile != NULL && decfile != NULL) {
		if (isroot) {
			cout << "Reading MPS file: " << mpsfile << endl;
			cout << "Reading DEC file: " << decfile << endl;
		}
		readMpsDec(env, mpsfile, decfile);
		//printModel(env);
	}

	if (paramfile != NULL) {
		if (isroot) cout << "Reading parameter files: " << paramfile << endl;
		readParamFile(env, paramfile);
	}

	if (isroot) cout << "Solving the problem with DW" << endl;
#ifdef DSP_HAS_MPI
	if (comm_size > 1)
		solveDwMpi(env, MPI_COMM_WORLD);
	else
#endif
		solveDw(env);

	if (isroot) cout << "Primal bound: " << getPrimalBound(env) << endl;

	if (isroot) cout << "Deleting DSP environment" << endl;
	freeEnv(env);
}

void readMpsDec(DspApiEnv* env, char* mpsfile, char* decfile) {
	// Read .mps file
	CoinMpsIO p;
	p.readMps(mpsfile);
	const CoinPackedMatrix* mps_matrix = p.getMatrixByRow();
	int ncols = p.getNumCols();
	const double* clbd = p.getColLower();
	const double* cubd = p.getColUpper();
	const double* obj = p.getObjCoefficients();
	vector<double> zeros(ncols, 0.0);
	string ctype;
	for (int j = 0; j < ncols; ++j)
		if (p.isContinuous(j))
			ctype.push_back('C');
		else
			ctype.push_back('I');

	// Read .dec file
	// For each block, the following vector stores the corresponding rows.
	vector<vector<string> > rows_in_blocks;
	parseDecFile(decfile, rows_in_blocks);

	// Assign block(s) to each process
	setBlockIds(env, rows_in_blocks.size() - 1);

	// For master block
	int blockid = 0;
	createBlockModel(env, p, mps_matrix, 0, rows_in_blocks[0], ctype.c_str(), obj);

	// For each sub-block
	for (unsigned i = 0; i < g_proc_idx_set.size(); ++i) {
		blockid = g_proc_idx_set[i] + 1; // blockid = subproblem index + 1, because master = 0
		//createBlockModel(env, p, mps_matrix, blockid, rows_in_blocks[blockid], ctype.c_str(), &zeros[0]);
		createBlockModel(env, p, mps_matrix, blockid, rows_in_blocks[blockid], ctype.c_str(), obj);
	}

	updateBlocks(env);
}

void parseDecFile(char* decfile, vector<vector<string> >& rows_in_blocks) {
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
	}
#if 0
	cout << "Number of blocks: " << rows_in_blocks.size() << endl;
	for (int i = 0; i < rows_in_blocks.size(); ++i) {
		cout << "Block " << i << endl;
		for (int j = 0; j < rows_in_blocks[i].size(); ++j)
			cout << rows_in_blocks[i][j] << endl;
	}
#endif
}

void createBlockModel(DspApiEnv* env, CoinMpsIO& p, const CoinPackedMatrix* mat, 
	int blockid, vector<string>& rows_in_block, const char* ctype, const double* obj) {
	vector<int> rowids;
	vector<double> rlbd, rubd;
	CoinPackedMatrix submat(false, 0, 0);

	for (unsigned j = 0; j < rows_in_block.size(); ++j) {
		for (int k = 0; k < p.getNumRows(); ++k) {
			if (strcmp(rows_in_block[j].c_str(), p.rowName(k)) == 0) {
				rowids.push_back(k);
				rlbd.push_back(p.getRowLower()[k]);
				rubd.push_back(p.getRowUpper()[k]);
				submat.appendRow(mat->getVectorSize(k), 
					mat->getIndices() + mat->getVectorFirst(k),
					mat->getElements() + mat->getVectorFirst(k));
				break;
			}
		}
	}
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
	loadBlockProblem(env, blockid, p.getNumCols(), rowids.size(), 
		submat.getNumElements(), submat.getVectorStarts(), submat.getIndices(), submat.getElements(), 
		p.getColLower(), p.getColUpper(), ctype, obj, &rlbd[0], &rubd[0]);
}

void setBlockIds(DspApiEnv* env, int nblocks, bool master_has_subblocks) {
	int comm_rank = 0;
	int comm_size = 1;
#ifdef DSP_HAS_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
#endif

	// empty block IDs
	g_proc_idx_set.clear();

	// DSP is further parallelized with comm_size > nblocks.
	int modrank = comm_rank % nblocks;
	// If we have more than one processor, do not assign a sub-block to the master.
	if (master_has_subblocks == false) {
		if (comm_size == 1) {
			// assign sub-blocks in round-robin fashion
			for (int s = modrank; s < nblocks; s += comm_size)
				g_proc_idx_set.push_back(s);
		} else {
			if (comm_rank > 0) {
				// exclude master
				comm_size--;
				modrank = (comm_rank-1) % nblocks;
	
				// assign sub-blocks in round-robin fashion
				for (int s = modrank; s < nblocks; s += comm_size) {
					g_proc_idx_set.push_back(s);
				}
			}
		}
	}

	// set the block ids to Dsp
	setIntPtrParam(env, "ARR_PROC_IDX", (int) g_proc_idx_set.size(), &g_proc_idx_set[0]);
}

