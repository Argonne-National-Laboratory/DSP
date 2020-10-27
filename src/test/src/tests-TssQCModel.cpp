// tests-BlkModel.cpp
#include "catch.hpp"

#include "DspCInterface.h"
#include "Model/TssModel.h"
#include "Model/QcModel.h"
#include "SolverInterface/DspOsiCpx.h"
#include "Solver/DualDecomp/DdDriverSerial.h"

void printCoreRows (DspApiEnv * env) {
    
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
    printCoreRows(env);
    
    /** test QcModel functions */
    bool isquadratic = true;
    setIsQuadratic(env, isquadratic);
    QcModel * qc = getQcModel(env);

    qc->readQuad(smpsfile, "../examples/quad/farmer");

    // qc->printQuadRows(0);
    // qc->printQuadRows(1);
    // qc->printQuadRows(2);

    // // x1 + 3x1^2 >= 1
    // int nqrows_0 = 1;
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
    // int nqrows_1 = 2;
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
    // int nqrows_2 = 2;
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

    // int nqrows[] = {nqrows_0, nqrows_1, nqrows_2};
    // qc->setQuadDimensions(3, nqrows);

    // qc->loadQuadraticConstrs(0, nqrows_0, linnzcnt_0, quadnzcnt_0, rhs_0, sense_0, linstart_0, linind_0, linval_0, quadstart_0, quadrow_0, quadcol_0, quadval_0);
    // qc->loadQuadraticConstrs(1, nqrows_1, linnzcnt_1, quadnzcnt_1, rhs_1, sense_1, linstart_1, linind_1, linval_1, quadstart_1, quadrow_1, quadcol_1, quadval_1);
    // qc->loadQuadraticConstrs(2, nqrows_2, linnzcnt_2, quadnzcnt_2, rhs_2, sense_2, linstart_2, linind_2, linval_2, quadstart_2, quadrow_2, quadcol_2, quadval_2);

    // qc->printQuadRows(0);
    // qc->printQuadRows(1);
    // qc->printQuadRows(2);
    
    /** set up solvers */
	// DSP_API_CHECK_MODEL();
	// freeSolver(env);

	// env->solver_ = new DdDriverSerial(env->model_, env->par_, env->message_);
	// DSP_RTN_CODE c = env->solver_->init();
    // REQUIRE(c == DSP_RTN_OK);
    
    for (int s = 0; s < getNumScenarios(env); s++)
    {
        QcRowDataScen * qcrowdata = qc->getQRowsCPXParams(s);

        // qc->printQuadRows(s);
        qc->printQuadRows(qcrowdata);

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
// 	// cpx->loadQuadraticConstrs(nqrows, linnzcnt, quadnzcnt, rhs, sense, linind, linval, quadrow, quadcol, quadval); 
//     cout << "reach 158" << endl;
//     qc->printQuadraticConstrs(s);

    }    
	

	// DSP_RTN_CHECK_THROW(dynamic_cast<DdDriverSerial*>(env->solver_)->run());
	// DSP_RTN_CHECK_THROW(env->solver_->finalize());


    freeEnv(env);
    REQUIRE(env == NULL);
}