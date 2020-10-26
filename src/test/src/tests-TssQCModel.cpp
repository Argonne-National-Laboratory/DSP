// tests-BlkModel.cpp
#include "catch.hpp"

#include "DspCInterface.h"
#include "Model/TssModel.h"
#include "Model/QcModel.h"
#include "SolverInterface/DspOsiCpx.h"
#include "Solver/DualDecomp/DdDriverSerial.h"

/** construct a map that maps variable names to their indices */
bool mapVarnameIndex(map<string, int> &map_varName_index, const char * corefilename) 
{

	ifstream corefile (corefilename);

	string name, item, rowname, varname;

	map_varName_index.clear();
    int nvars = 0;
	vector<string> rowNames;

	map<string, int>::iterator map_it;
    vector<string>::iterator vec_it;

    bool foundColumns = false;

    if (corefile.is_open()) {
        
	    while (corefile >> item) {
		
            if (item.find("ROWS") != string::npos) {
			
                while (1) {
			        corefile >> item >> name;

                    if (item.find("COLUMNS") != string::npos) {
                        foundColumns = true;
                        map_varName_index[name] = nvars;
                        nvars++;
                        break;
                    }

                    if (name == "ENDATA" || item == "ENDATA") {
                        cout << "Encountered ENDATA before reading Columns" << endl;
                        break;
                    }

			        rowNames.push_back(name);
                }
		
                if (foundColumns) {
                
                    while(1) {
				
	        			corefile >> name >> item;
				
			        	vec_it = find(rowNames.begin(), rowNames.end(), name);
				        if (vec_it == rowNames.end()) 
				        {
					        /** var term */
					        corefile >> item;

					        if (name.find("RHS") != string::npos)
						        break;
					
					        map_it = map_varName_index.find(name);
  					        if (map_it == map_varName_index.end())
    					    {
    					        map_varName_index[name] = nvars;
                                nvars++;
					        }
				        } else 
				        {
					        /* constraint term
					        * does not do anything */
				        }
			        }
		        } else {
                    cout << "Columns data should follow Rows data" << endl;
                    break;
                }
            }

		if (item.find("ENDATA") != string::npos) 
			break;
        
	    }
    corefile.close();
    } else {
        cout << "Unable to open file";
        return false;
    }
                             
	return true;
}

void printCoreRows (DspApiEnv * env, map<string, int> &map_varName_index) {
    TssModel * tss = dynamic_cast<TssModel*>(env->model_);
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
    if (mapVarnameIndex(map_varName_index, "../examples/smps/farmer.cor")) 
    {
        cout << "variable, index" << endl;
        for (auto &v : map_varName_index)
            cout << v.first << ", " << v.second << endl;
    } else cout << "fail to get a list of varNames" << endl;
    
}

TEST_CASE("Class TssQCModel") {
    
    DspApiEnv* env = NULL;

    /** createEnv function */
    env = createEnv();
    REQUIRE(env != NULL);
    REQUIRE(env->model_ == NULL);
    REQUIRE(env->solver_ == NULL);

    /** read SMPS file */        
    char smpsfile[80];
    sprintf(smpsfile, "../examples/smps/farmer");
    puts(smpsfile);
            
    int ret = readSmps(env, smpsfile);
    REQUIRE(ret == 0);

    /** test variable indices by printing constraints */
    map<string, int> map_varName_index;
    printCoreRows(env, map_varName_index);
    REQUIRE(map_varName_index.size() == getNumCols(env, 0) + getNumCols(env, 1));

    /** test QcModel functions */
    QcModel * qc;
	try
	{
		qc = dynamic_cast<QcModel *>(env->model_);
	} catch (const std::bad_cast& e)
	{
		printf("Error: Model claims to be QC when it is not");
	}

    qc->readQuad(map_varName_index, "../examples/quad/farmer.txt");

    qc->printQuadraticConstrs(0);
    qc->printQuadraticConstrs(1);
    qc->printQuadraticConstrs(2);

    // // x1 + 3x1^2 >= 1
    // int nqconstrs_0 = 1;
    // int linnzcnt_0[] = {1};
    // int quadnzcnt_0[] = {1};
    // double rhs_0[] = {1};
    // int sense_0[] = {'G'};
    // int linstart_0[] = {0};
    // int quadstart_0[] = {0};
    
    // int linind_0[] = {1};
    // double linval_0[] = {1};
    // int quadrow_0[] = {1};
    // int quadcol_0[] = {1};
    // double quadval_0[] = {3};

    // // x2 + x2^2 >= 1
    // // x1x2 >= 3
    // /** qconstrs 갯수 만큼 */
    // int nqconstrs_1 = 2;
    // int linnzcnt_1[] = {1, 0};
    // int quadnzcnt_1[] = {1, 1};
    // double rhs_1[] = {1, 3};
    // int sense_1[] = {'G', 'G'};
    // int linstart_1[] = {0,0}; 
    // int quadstart_1[] = {0,1}; 
    // // *** /
    // int linind_1[] = {2};
    // double linval_1[] = {1};
    // int quadrow_1[] = {2, 1};
    // int quadcol_1[] = {2, 2};
    // double quadval_1[] = {1, 1};

    // // 5x1 + 4 x2 + 2x1x2 + 3x1^2 + 2x^2 >= 2
    // // 2x1 + 2x1x2 + 2x^2 >= 5
    // int nqconstrs_2 = 2;
    // /** qconstrs 갯수 만큼 */
    // int linnzcnt_2[] = {2,1};
    // int quadnzcnt_2[] = {3,2};
    // double rhs_2[] = {2,5};
    // int sense_2[] = {'G','G'};
    // int linstart_2[] = {0,2};
    // int quadstart_2[] = {0,3};
    // // *** /
    // int linind_2[] = {1, 2,1,2};
    // double linval_2[] = {5, 4,2,2};
    // int quadrow_2[] = {1, 1, 2,1,2};
    // int quadcol_2[] = {2, 1, 2,1,2};
    // double quadval_2[] = {2, 3, 2,2,2};

    // int nqconstrs[] = {nqconstrs_0, nqconstrs_1, nqconstrs_2};
    // qc->setQuadDimensions(3, nqconstrs);

    // qc->loadQuadraticConstrs(0, nqconstrs_0, linnzcnt_0, quadnzcnt_0, rhs_0, sense_0, linstart_0, linind_0, linval_0, quadstart_0, quadrow_0, quadcol_0, quadval_0);
    // qc->loadQuadraticConstrs(1, nqconstrs_1, linnzcnt_1, quadnzcnt_1, rhs_1, sense_1, linstart_1, linind_1, linval_1, quadstart_1, quadrow_1, quadcol_1, quadval_1);
    // qc->loadQuadraticConstrs(2, nqconstrs_2, linnzcnt_2, quadnzcnt_2, rhs_2, sense_2, linstart_2, linind_2, linval_2, quadstart_2, quadrow_2, quadcol_2, quadval_2);

    // qc->printQuadraticConstrs(0);
    // qc->printQuadraticConstrs(1);
    // qc->printQuadraticConstrs(2);
    
    // /** set up solvers */
	// DSP_API_CHECK_MODEL();
	// freeSolver(env);

	// env->solver_ = new DdDriverSerial(env->model_, env->par_, env->message_);
	// DSP_RTN_CODE c = env->solver_->init();
    // REQUIRE(c == DSP_RTN_OK);
    
    // for (int s = 0; s < getNumScenarios(env); s++)
    // {
    //     /** quadratic constraints data */
	//     int nqconstrs = qc->getNumQConstrs(s);
	//     int * linnzcnt = qc->getLinearNonZeroCounts(s);
	//     int * quadnzcnt = qc->getQuadNonZeroCounts(s); 
	//     double * rhs = qc->getRhs(s); 
	//     int * sense = qc->getSense(s); 
	//     const int ** linind = qc->getLinearIndices(s); 
	//     const int ** quadrow = qc->getQuadIndices1st(s); 
	//     const int ** quadcol = qc->getQuadIndices2nd(s); 
	//     const double ** linval = qc->getLinearVals(s); 
	//     const double ** quadval = qc->getQuadraticVals(s); 

    //     qc->printQuadraticConstrs(s);

//         for (int i = 0; i < nqconstrs; i++) 
// 	{
// 		cout << "Scen " << s << "th " << i << "th quad constr: ";

// 		for (int lt = 0; lt < linnzcnt[i]; lt++)
// 		{
// 			cout << linval[i][lt] << " x" << linind[i][lt] << " + ";
// 		}
// 		for (int qt = 0; qt < quadnzcnt[i]-1; qt++)
// 		{
// 			cout << quadval[i][qt] << " x" << quadrow[i][qt] << " x" << quadcol[i][qt] << " + ";
// 		}
// 		cout << quadval[i][quadnzcnt[i]-1] << " x" << quadrow[i][quadnzcnt[i]-1] << " x" << quadcol[i][quadnzcnt[i]-1];
// 		if (sense[i] == 'L')
// 			cout << " <= " << rhs[i] << endl;
// 		else 
// 			cout << " >= " << rhs[i] << endl;

// 	}
//     cout << "reach 154" << endl;
//     DdWorkerLB* worker_ = dynamic_cast<DdMWSerial*> (env->solver_)

//     /** retrieve DdWorker pointers */
// 	for (unsigned i = 0; i < worker_.size(); ++i)
// 	{
// 		switch(worker_[i]->getType())
// 		{
// 		case DdWorker::LB:
// 			workerlb = dynamic_cast<DdWorkerLB*>(worker_[i]);
// 			break;
// 		case DdWorker::CGBd:
// #ifdef DSP_HAS_SCIP
// 			workercg = dynamic_cast<DdWorkerCGBd*>(worker_[i]);
// #endif
// 			break;
// 		case DdWorker::UB:
// 			workerub = dynamic_cast<DdWorkerUB*>(worker_[i]);
// 			break;
// 		default:
// 			message_->print(0, "Unknown worker type (%d).\n", worker_[i]->getType());
// 			break;
// 		}
// 	}

//     // DspOsiCpx * cpx = dynamic_cast<DspOsiCpx*> (dynamic_cast<DdDriverSerial*> (env->solver_));
//     cout << "reach 156" << endl;
// 	// cpx->loadQuadraticConstrs(nqconstrs, linnzcnt, quadnzcnt, rhs, sense, linind, linval, quadrow, quadcol, quadval); 
//     cout << "reach 158" << endl;
//     qc->printQuadraticConstrs(s);

    //     FREE_ARRAY_PTR(rhs);
    //     FREE_ARRAY_PTR(sense);
    //     FREE_ARRAY_PTR(linnzcnt);
    //     FREE_ARRAY_PTR(quadnzcnt);  
    //     // FREE_2D_ARRAY_PTR(nqconstrs,linind);
    //     // FREE_2D_ARRAY_PTR(nqconstrs,linval);
	//     // FREE_2D_ARRAY_PTR(nqconstrs,quadrow);
	//     // FREE_2D_ARRAY_PTR(nqconstrs,quadcol);
	//     // FREE_2D_ARRAY_PTR(nqconstrs,quadval);
    // }    
	

	// DSP_RTN_CHECK_THROW(dynamic_cast<DdDriverSerial*>(env->solver_)->run());
	// DSP_RTN_CHECK_THROW(env->solver_->finalize());


    freeEnv(env);
    REQUIRE(env == NULL);
}