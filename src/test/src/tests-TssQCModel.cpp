// tests-BlkModel.cpp
#include "catch.hpp"

#include "DspCInterface.h"
#include "Model/TssModel.h"
#include "Model/DecTssQcModel.h"

int runDsp(char* algotype, char* smpsfile, char* mpsfile, char* decfile, char* solnfile, char* paramfile, char* testvalue, char* quadfile);
int readMpsDec(DspApiEnv* env, char* mpsfile, char* decfile);
int parseDecFile(char* decfile, vector<vector<string> >& rows_in_blocks);
void setBlockIds(DspApiEnv* env, int nsubprobs, bool master_has_subblocks);
void createBlockModel(DspApiEnv* env, CoinMpsIO& p, const CoinPackedMatrix* mat, 
	int blockid, vector<string>& rows_in_block, map<string,int>& rowname2index, 
	const char* ctype, const double* obj);
const double test_tolerance = 1.0e-2;
void printCoreRows (DspApiEnv * env);

void copyToCharArray (string &stringToCopy, char* &arrayToWrite) {
    
	if (arrayToWrite != NULL) {
		free(arrayToWrite);
		arrayToWrite = NULL;
	}
    arrayToWrite = new char [stringToCopy.length()];
    strcpy(arrayToWrite, stringToCopy.c_str());
}

TEST_CASE("Test runDSP") {
	
	bool isroot = true;
/* test with MPI turned off */
#undef DSP_HAS_MPI

    char* algotype = NULL;
    char* smpsfile = NULL;
    char* mpsfile = NULL;
    char* decfile = NULL;
    char* solnfile = NULL;
    char* paramfile = NULL;
    char* testvalue = NULL;
    char* quadfile = NULL;

	string algotype_s;
	string smpsfile_s;
	string mpsfile_s;
	string decfile_s;
	string solnfile_s;
	string paramfile_s;
	string testvalue_s;
	string quadfile_s;

#ifdef DSP_HAS_SCIP

    paramfile_s = "../test/params_scip.txt";
    copyToCharArray(paramfile_s, paramfile);

    SECTION ("Solving a MPS-DEC Instance") {

        mpsfile_s = "../examples/mps-dec/noswot.mps";
        decfile_s = "../examples/mps-dec/noswot.dec";
		testvalue_s = "-42";
        copyToCharArray(mpsfile_s, mpsfile);
        copyToCharArray(decfile_s, decfile);
		copyToCharArray(testvalue_s, testvalue);

        SECTION("with DE") {
            algotype_s = "de";
            copyToCharArray(algotype_s, algotype);
            
            REQUIRE(runDsp(algotype, smpsfile, mpsfile, decfile, solnfile, paramfile, testvalue, quadfile) == 0);

			delete[] algotype;
        }

        SECTION("with DD") {
            algotype_s = "dd";
            copyToCharArray(algotype_s, algotype);
            
            REQUIRE(runDsp(algotype, smpsfile, mpsfile, decfile, solnfile, paramfile, testvalue, quadfile) == 0);

			delete[] algotype;
        }

        // SECTION("DW") {
        //     algotype_s = "dw";
        //     copyToCharArray(algotype_s, algotype);
                
        //     REQUIRE(runDsp(algotype, smpsfile, mpsfile, decfile, solnfile, paramfile, testvalue, quadfile) == 0);
        // }
		
		delete[] mpsfile;
		delete[] decfile;
		delete[] testvalue;
    }

    SECTION ("Solving a SMPS Instance") {

        smpsfile_s = "../examples/smps/farmer";
		copyToCharArray(smpsfile_s, smpsfile);

		SECTION("with DE") {
            algotype_s = "de";
			testvalue_s = "-108389.9994043";
            copyToCharArray(algotype_s, algotype);
			copyToCharArray(testvalue_s, testvalue);
            
            REQUIRE(runDsp(algotype, smpsfile, mpsfile, decfile, solnfile, paramfile, testvalue, quadfile) == 0);

			delete[] algotype;
			delete[] testvalue;
        }

        SECTION("with DD") {
            algotype_s = "dd";
			testvalue_s = "-108389.9994043";
            copyToCharArray(algotype_s, algotype);
			copyToCharArray(testvalue_s, testvalue);
            
            REQUIRE(runDsp(algotype, smpsfile, mpsfile, decfile, solnfile, paramfile, testvalue, quadfile) == 0);

			delete[] algotype;
			delete[] testvalue;
        }

        SECTION("with Quadratic Constraints") {
			
			testvalue_s = "-105093";
			copyToCharArray(testvalue_s, testvalue);

            quadfile_s = "../examples/quad/farmer";
            copyToCharArray(quadfile_s, quadfile);
			puts(quadfile);

            SECTION("with DE") {
				algotype_s = "de";
				copyToCharArray(algotype_s, algotype);
				
				REQUIRE(runDsp(algotype, smpsfile, mpsfile, decfile, solnfile, paramfile, testvalue, quadfile) == 0);

				delete[] algotype;
			}

			SECTION("with DD") {
				algotype_s = "dd";
				copyToCharArray(algotype_s, algotype);
				
				REQUIRE(runDsp(algotype, smpsfile, mpsfile, decfile, solnfile, paramfile, testvalue, quadfile) == 0);

				delete[] algotype;
			}
			
			delete[] quadfile;
			delete[] testvalue;
        }
		delete[] smpsfile;
    }

	delete[] paramfile;
#endif

#ifdef DSP_HAS_CPX

    paramfile_s = "../test/params_cpx.txt";
    copyToCharArray(paramfile_s, paramfile);

    SECTION ("Solving a MPS-DEC Instance") {

        mpsfile_s = "../examples/mps-dec/noswot.mps";
        decfile_s = "../examples/mps-dec/noswot.dec";
		testvalue_s = "-42";
        copyToCharArray(mpsfile_s, mpsfile);
        copyToCharArray(decfile_s, decfile);
		copyToCharArray(testvalue_s, testvalue);

        SECTION("with DE") {
            algotype_s = "de";
            copyToCharArray(algotype_s, algotype);
            
            REQUIRE(runDsp(algotype, smpsfile, mpsfile, decfile, solnfile, paramfile, testvalue, quadfile) == 0);

			delete[] algotype;
        }

        SECTION("with DD") {
            algotype_s = "dd";
            copyToCharArray(algotype_s, algotype);
            
            REQUIRE(runDsp(algotype, smpsfile, mpsfile, decfile, solnfile, paramfile, testvalue, quadfile) == 0);

			delete[] algotype;
        }

        // SECTION("DW") {
        //     algotype_s = "dw";
        //     copyToCharArray(algotype_s, algotype);
                
        //     REQUIRE(runDsp(algotype, smpsfile, mpsfile, decfile, solnfile, paramfile, testvalue, quadfile) == 0);
        // }
		delete[] mpsfile;
		delete[] decfile;
		delete[] testvalue;
    }

    SECTION ("Solving a SMPS Instance") {

        smpsfile_s = "../examples/smps/farmer";
		copyToCharArray(smpsfile_s, smpsfile);
		
		SECTION("with DE") {
            algotype_s = "de";
			testvalue_s = "-108389.9994043";
            copyToCharArray(algotype_s, algotype);
			copyToCharArray(testvalue_s, testvalue);
            
            REQUIRE(runDsp(algotype, smpsfile, mpsfile, decfile, solnfile, paramfile, testvalue, quadfile) == 0);

			delete[] algotype;
			delete[] testvalue;
        }

        SECTION("with DD") {
            algotype_s = "dd";
			testvalue_s = "-108389.9994043";
            copyToCharArray(algotype_s, algotype);
			copyToCharArray(testvalue_s, testvalue);
            
            REQUIRE(runDsp(algotype, smpsfile, mpsfile, decfile, solnfile, paramfile, testvalue, quadfile) == 0);

			delete[] algotype;
			delete[] testvalue;
        }

        SECTION("with Quadratic Constraints") {
			
			SECTION("farmer1") {
				testvalue_s = "-105093";
				copyToCharArray(testvalue_s, testvalue);

				quadfile_s = "../examples/quad/farmer";
				copyToCharArray(quadfile_s, quadfile);
				puts(quadfile);

				SECTION("with DE") {
					algotype_s = "de";
					copyToCharArray(algotype_s, algotype);
					
					REQUIRE(runDsp(algotype, smpsfile, mpsfile, decfile, solnfile, paramfile, testvalue, quadfile) == 0);

					delete[] algotype;
				}

				SECTION("with DD") {
					algotype_s = "dd";
					copyToCharArray(algotype_s, algotype);
					
					REQUIRE(runDsp(algotype, smpsfile, mpsfile, decfile, solnfile, paramfile, testvalue, quadfile) == 0);

					delete[] algotype;
				}

				delete[] testvalue;
				delete[] quadfile;
			}

			SECTION("farmer2") {
				testvalue_s = "-44147.4";
				copyToCharArray(testvalue_s, testvalue);

				quadfile_s = "../examples/quad/farmer2";
				copyToCharArray(quadfile_s, quadfile);
				puts(quadfile);

				SECTION("with DE") {
					algotype_s = "de";
					copyToCharArray(algotype_s, algotype);
					
					REQUIRE(runDsp(algotype, smpsfile, mpsfile, decfile, solnfile, paramfile, testvalue, quadfile) == 0);

					delete[] algotype;
				}

				SECTION("with DD") {
					algotype_s = "dd";
					copyToCharArray(algotype_s, algotype);
					
					REQUIRE(runDsp(algotype, smpsfile, mpsfile, decfile, solnfile, paramfile, testvalue, quadfile) == 0);

					delete[] algotype;
				}

				delete[] testvalue;
				delete[] quadfile;
			}
			
        }
		delete[] smpsfile;
    }
	delete[] paramfile;
#endif

#ifdef DSP_HAS_MPI
		MPI_Finalize();
#endif
}


int runDsp(char* algotype, char* smpsfile, char* mpsfile, char* decfile, char* solnfile, char* paramfile, char* testvalue, char* quadfile) {

	int ret = 0;
	bool isroot = true;
	bool issolved = true;
	bool isstochastic = smpsfile != NULL ? true : false;
	bool isquadratic = quadfile != NULL ? true : false;
#ifdef DSP_HAS_MPI
	int comm_rank, comm_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
	isroot = comm_rank == 0 ? true : false;
#endif

	if (isroot) cout << "Creating DSP environment\n";
	DspApiEnv* env = createEnv();

	/* create model */
	ret = createModel(env, isstochastic, isquadratic);
	if (ret != 0) return ret;

	// Read problem instance from file(s)
	if (smpsfile != NULL) 
	{
		if (isroot) cout << "Reading SMPS files: " << smpsfile << endl;
		ret = readSmps(env, smpsfile);
		if (ret != 0) return ret;
		if (isroot) 
		{
			cout << "First stage: " << getNumRows(env,0) << " rows, " << getNumCols(env,0) << " cols, " << getNumIntegers(env,0) << " integers" << endl;
			cout << "Second stage: " << getNumRows(env,1) << " rows, " << getNumCols(env,1) << " cols, " << getNumIntegers(env,1) << " integers" << endl;
			cout << "Number of scenarios: " << getNumSubproblems(env) << endl;
		}
		setBlockIds(env, getNumSubproblems(env), true);
	} else if (mpsfile != NULL && decfile != NULL) 
	{
		if (isroot) 
		{
			cout << "Reading MPS file: " << mpsfile << endl;
			cout << "Reading DEC file: " << decfile << endl;
		}
		ret = readMpsDec(env, mpsfile, decfile);
		if (ret != 0) return ret;
		isstochastic = false;
	}
	
		if (isroot) 
	if (paramfile != NULL) {
		if (isroot) cout << "Reading parameter files: " << paramfile << endl;
		readParamFile(env, paramfile);
	}

	if (quadfile != NULL) 
	{
		if (string(algotype) != "de" && string(algotype) != "dd") {
			cout << "Current version only support deterministic or dual decomposition solvers for quadratic constrained problem" << endl;
			return 1;
		}

		#ifndef DSP_HAS_CPX
			cout << "Current version only support CPLEX for solving quadratic constrained problem" << endl;
			return 1;
		#endif
		/** force to use CPLEX if available */
		env->par_->setIntParam("DD/MASTER/SOLVER", OsiCpx);
		env->par_->setIntParam("DD/SUB/SOLVER", OsiCpx);
		env->par_->setIntParam("DE/SOLVER", OsiCpx);

		if (isroot) cout << "Reading Quad files: " << quadfile << endl;
		ret = readQuad(env, smpsfile, quadfile);
		if (ret != 0) return ret;
		if (isroot) 
		{
			for (int s = 0; s < getNumScenarios(env); s++)
				cout << "Second stage: " << getNumQRows(env,s) << " quadratic rows in scenario " << s << endl;
		}
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
					if ((primobj - val) / (fabs(val) + 1.e-10) > test_tolerance || (val - dualobj) / (fabs(val) + 1.e-10) > test_tolerance)
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

	// cout << "Creating hash table for matrix rows ... ";
	map<string, int> rowname2index;
	for (int i = 0; i < p.getNumRows(); ++i)
		rowname2index[string(p.rowName(i))] = i;
	// cout << "done!" << endl;

	// Read .dec file
	// For each block, the following vector stores the corresponding rows.
	// cout << "Parsing .dec file ... ";
	vector<vector<string> > rows_in_blocks;
	ret = parseDecFile(decfile, rows_in_blocks);
	if (ret != 0) return ret;

	// cout <<  "Assign block(s) to each process" << endl;
	setBlockIds(env, rows_in_blocks.size() - 1, true);

	vector<int> proc_idx_set;
	for (int s = 0; s < env->par_->getIntPtrParamSize("ARR_PROC_IDX"); ++s) {
		proc_idx_set.push_back(env->par_->getIntPtrParam("ARR_PROC_IDX")[s]);
	}

	// cout << "create master block" << endl;
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

void printCoreRows (DspApiEnv * env) {
    
    TssModel * tss = getTssModel(env);
	cout << "number of stages: " << tss->getNumStages() << endl;
	cout << "row data: " << endl;
	    
    int i, j, k, s;
	for (s = 0; s < tss->getNumStages(); ++s)
	{
        cout << "stage " << s << ": " << endl; 
		
        double * rlbd_core = new double [tss->getNumRows(s)];
        double * rubd_core = new double [tss->getNumRows(s)];
        tss->copyCoreRowLower(rlbd_core, s);
        tss->copyCoreRowUpper(rubd_core, s);

        for (j = 0; j < tss->getNumRows(s); ++j)
		{
            if (s == 0) i = j;
            else i = j + tss->getNumRows(s-1);

            const CoinPackedVector * row = tss->getRowCore(i);

            cout << "constr " << j << " : " << rlbd_core[j] << " <= ";

            int nTerms = row->getNumElements();

            for (k = 0; k < nTerms-1; k++) {
				cout << row->getElements()[k] << " * " << row->getIndices()[k] << " + ";
			}
		    cout << row->getElements()[nTerms-1] << " * " << row->getIndices()[nTerms-1] << " <= " << rubd_core[j];
			cout << endl;
		}
	}
}