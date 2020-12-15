#include "catch.hpp"

#include "DspCInterface.h"
#include "Model/TssModel.h"

/** data for farmer example with quadratic objective */

int NS = 3;                                                            // number of scenarios
double probability[3] = {(double)1 / 3, (double)1 / 3, (double)1 / 3}; // probability

/** first-stage constraint matrix */
CoinBigIndex start1[2] = {0, 3};
const int index1[3] = {0, 1, 2};
const double value1[3] = {1, 1, 1};

/** quadratic objective information for firtst stage:
*  i.e., C0^2+2*C0*C1+5*C1^2
*/
const int qrowindex[3] = {0, 0, 1};
const int qcolindex[3] = {0, 1, 1};
const double qvalue[3] = {1, 2, 5};
const CoinBigIndex numq = 3;
//int qindice1[3]={0, 1, 1};
int qstart1[4]={0, 2, 3, 3};
//first-stage variable lower bound, upper bound, and type
const double xclbd[3] = {0, 0, 0};
const double xcubd[3] = {500, 500, 500};
const char xctype[3] = {'C', 'C', 'C'};

//first-stage objective
const double Cost[3] = {150, 230, 260}; // cost of planting crops

//first-stage constraint upper and lower limit
const double rubd1[1] = {500}; // budget capacity
const double rlbd1[1] = {0};   // budget capacity

/** second-stage constraint matrix
 *  the index starts from the first-stage variables,
 *  and the maximum is less than 
 *  # first-stage variables + # of single scenario second-stage variables
 */

/** second-stage objective coefficient: only containes second stage variables
 *   selling price 238 210
 *  purchase price -170, -150, -36, -10
 */
const double profit[6] = {238, 210, -170, -150, -36, -10};

/** second-stage constraint matrix */
const CoinBigIndex start[4] = {0, 3, 6, 9};
const int index2[9] = {0, 3, 5, 1, 4, 6, 2, 7, 8};
const double value[3][9] = {
    {3.0, 1, -1, 3.6, 1, -1, 24.0, -1, -1},
    {2.5, 1, -1, 3.0, 1, -1, 20.0, -1, -1},
    {2.0, 1, -1, 2.4, 1, -1, 16.0, -1, -1}};

/** lower and upper limits of second-stage variables */
const double yclbd[6] = {0, 0, 0, 0, 0, 0};
const double ycubd[6] = {INFINITY, INFINITY, INFINITY, INFINITY, 6000, INFINITY};

//second stage quadratic objective
    const int qrowindex2[3][3] = 
    {   {3, 3, 4},
        {3, 4, 5},
        {3, 3, 5}};
    const int qcolindex2[3][3] = {
        {3, 4, 4},
        {3, 4, 5},
        {3, 5, 5}
        };
    const double qvalue2[3][3] = {
        {1, 2, 5},
        {1, 1, 2},
        {2, 0.5, 4}
        };
int qstart2[3][10]={
    {0, 0, 0, 0, 2, 3, 3, 3, 3, 3},
    {0, 0, 0, 0, 1, 2, 3, 3, 3, 3},
    {0, 0, 0, 0, 2, 2, 3, 3, 3, 3}
};
// const int qrowindex2[3] = {3, 3, 4};
// const int qcolindex2[3] = {3, 4, 4};
// const double qvalue2[3] = {1, 2, 5};
const CoinBigIndex numq2 = 3;
//variable type of second-stage variables
const char yctype[6] = {'C', 'C', 'C', 'C', 'C', 'C'};

//lower and upper limits of constraints: minimum crop requirement
const double rlbd[3] = {200, 240, 0};
const double rubd[3] = {INFINITY, INFINITY, INFINITY};

TEST_CASE("Check the parameters loaded by functions")
{
    DspApiEnv *env = createEnv();
    //createModel(env, 1, 0);
    /** set number of scenarios */

    setNumberOfScenarios(env, NS);

    /** set problem dimension */
    setDimensions(env, 3, 1, 6, 3);

    SECTION("load first stage")
    {
        SECTION("with quadratic objectives in the 1st stage")
        {
            loadQuadraticFirstStage(env, start1, index1, value1, xclbd, xcubd, xctype, Cost, qrowindex, qcolindex, qvalue, numq, rubd1, rubd1);
            TssModel *model = getTssModel(env);
            /** two-stage problem 
              * check number of stages;
              * check number of scenarios;
            */
            REQUIRE(model->getNumStages()==2);
            REQUIRE(model->getNumScenarios()==NS);

            /** number of columns in the first stage problem */
            int numcol1=model->getNumCols(0);
            REQUIRE(numcol1==3);

            /** check first stage core[0] */
            const double *objcore=new double[numcol1];
            objcore=model->getObjCore(0);
            for (int i=0; i<numcol1; i++){
                REQUIRE(objcore[i]==Cost[i]);
            }
            /** check quadratic objectives in the first stage */
            const CoinPackedMatrix *qobj1;
            qobj1=model->getQuadraticObjCore(0);
            REQUIRE(qobj1->getMajorDim()==numcol1);
            for (int i=0; i<qobj1->getNumElements(); i++){
                REQUIRE(qobj1->getElements()[i]==qvalue[i]);
                REQUIRE(qobj1->getIndices()[i]==qcolindex[i]);
            }
            for (int i=0; i<qobj1->getSizeVectorStarts(); i++){
                REQUIRE(qobj1->getVectorStarts()[i]==qstart1[i]);
            }

            SECTION("load second stage with stochastic T and W")
            {
                SECTION("with quadratic objectives in the 2nd stage")
                {
                    for (int i = 0; i < NS; i++)
                    {
                        /** the objective function is 
                            150 C0 + 230 C1 + 260 C2 + 79.33333333333333 C3 + 70 C4
                            - 56.66666666666666 C5 - 50 C6 - 12 C7 - 3.333333333333333 C8
                            + 79.33333333333333 C9 + 70 C10 - 56.66666666666666 C11 - 50 C12
                            - 12 C13 - 3.333333333333333 C14 + 79.33333333333333 C15 + 70 C16
                            - 56.66666666666666 C17 - 50 C18 - 12 C19 - 3.333333333333333 C20 + [
                            2 C0 ^2 + 4 C0 * C1 + 10 C1 ^2 + 0.6666666666666666 C3 ^2
                            + 1.333333333333333 C3 * C4 + 3.333333333333333 C4 ^2
                            + 0.6666666666666666 C9 ^2 + 1.333333333333333 C9 * C10
                            + 3.333333333333333 C10 ^2 + 0.6666666666666666 C15 ^2
                            + 1.333333333333333 C15 * C16 + 3.333333333333333 C16 ^2 ] / 2 
                        */
                        loadQuadraticSecondStage(env, i, probability[i], start, index2, value[i], yclbd, ycubd, yctype, profit, qrowindex2[0], qcolindex2[0], qvalue2[0], numq2, rlbd, rubd);
                    }
                    TssModel *model = getTssModel(env);
                    
                    /** check linear objectives for each scenario in the second stage */
                    const double * elements1;
                    for (int i=0; i<model->getNumScenarios(); i++){
                        elements1=new double [model->getObjScenario(i)->getNumElements()];
                        elements1=model->getObjScenario(i)->getElements();
                        for (int j=0; j<model->getObjScenario(i)->getNumElements(); j++){
                            REQUIRE(elements1[i]==profit[i]);
                        }
                    }
                    
                    
                    /** check quadratic objectives in the second stage */
                    const CoinPackedMatrix *qobj2;
                    
                    for (int i=0; i<model->getNumScenarios(); i++){
                        qobj2=model->getQuadraticObjScenario(i);
                        for (int j=0; j<qobj2->getNumElements(); j++){
                            REQUIRE(qobj2->getElements()[j]==qvalue2[0][j%3]);
                            REQUIRE(qobj2->getIndices()[j]==qcolindex2[0][j%3]);
                        }
        
                        for (int j=0; j<qobj2->getSizeVectorStarts(); j++){
                            REQUIRE(qobj2->getVectorStarts()[j]==qstart2[0][j]);
                        }
                    }
                }
            }
        }
    }
}

TEST_CASE("Extensive Formulation for MIQP")
{
    DspApiEnv *env = createEnv();
    //createModel(env, 1, 0);
    /** set number of scenarios */
    setNumberOfScenarios(env, NS);

    /** set problem dimension */
    setDimensions(env, 3, 1, 6, 3);

    SECTION("load first stage")
    {
        SECTION("with quadratic objectives in the 1st stage")
        {
            loadQuadraticFirstStage(env, start1, index1, value1, xclbd, xcubd, xctype, Cost, qrowindex, qcolindex, qvalue, numq, rubd1, rubd1);
            TssModel *model = getTssModel(env);
            /** two-stage problem 
              * check number of stages;
              * check number of scenarios;
            */
            REQUIRE(model->getNumStages()==2);
            REQUIRE(model->getNumScenarios()==NS);

            /** number of columns in the first stage problem */
            int numcol1=model->getNumCols(0);
            REQUIRE(numcol1==3);

            /** check linear objectives in each stage */
            /** check first stage core[0] */
            const double *objcore=new double[numcol1];
            objcore=model->getObjCore(0);
            for (int i=0; i<numcol1; i++){
                REQUIRE(objcore[i]==Cost[i]);
            }
            
            SECTION("load second stage with stochastic T and W")
            {
                SECTION("with quadratic objectives in the 2nd stage")
                {
                    for (int i = 0; i < NS; i++)
                    {
                        /** the objective function is 
                            150 C0 + 230 C1 + 260 C2 + 79.33333333333333 C3 + 70 C4
                            - 56.66666666666666 C5 - 50 C6 - 12 C7 - 3.333333333333333 C8
                            + 79.33333333333333 C9 + 70 C10 - 56.66666666666666 C11 - 50 C12
                            - 12 C13 - 3.333333333333333 C14 + 79.33333333333333 C15 + 70 C16
                            - 56.66666666666666 C17 - 50 C18 - 12 C19 - 3.333333333333333 C20 + [
                            2 C0 ^2 + 4 C0 * C1 + 10 C1 ^2 + 0.6666666666666666 C3 ^2
                            + 1.333333333333333 C3 * C4 + 3.333333333333333 C4 ^2
                            + 0.6666666666666666 C9 ^2 + 1.333333333333333 C9 * C10
                            + 3.333333333333333 C10 ^2 + 0.6666666666666666 C15 ^2
                            + 1.333333333333333 C15 * C16 + 3.333333333333333 C16 ^2 ] / 2 
                        */
                        loadQuadraticSecondStage(env, i, probability[i], start, index2, value[i], yclbd, ycubd, yctype, profit, qrowindex2[0], qcolindex2[0], qvalue2[0], numq2, rlbd, rubd);
                    }

                    

                    /** check linear objectives for each scenario in the second stage */
                    const double * elements1;
                    for (int i=0; i<model->getNumScenarios(); i++){
                        elements1=new double [model->getObjScenario(i)->getNumElements()];
                        elements1=model->getObjScenario(i)->getElements();
                        for (int j=0; j<model->getObjScenario(i)->getNumElements(); j++){
                            REQUIRE(elements1[i]==profit[i]);
                        }
                    }
                    
                    /** check quadratic objectives in the first stage */
                    const CoinPackedMatrix *qobj1;
                    qobj1=model->getQuadraticObjCore(0);
                    REQUIRE(qobj1->getMajorDim()==numcol1);
                    for (int i=0; i<qobj1->getNumElements(); i++){
                        REQUIRE(qobj1->getElements()[i]==qvalue[i]);
                        REQUIRE(qobj1->getIndices()[i]==qcolindex[i]);
                    }
                    for (int i=0; i<qobj1->getSizeVectorStarts(); i++){
                        REQUIRE(qobj1->getVectorStarts()[i]==qstart1[i]);
                    }

                    /** check quadratic objectives in the second stage */
                    const CoinPackedMatrix *qobj2;
                    
                    for (int i=0; i<model->getNumScenarios(); i++){
                        qobj2=model->getQuadraticObjScenario(i);
                        for (int j=0; j<qobj2->getNumElements(); j++){
                            REQUIRE(qobj2->getElements()[j]==qvalue2[0][j%3]);
                            REQUIRE(qobj2->getIndices()[j]==qcolindex2[0][j%3]);
                        }
        
                        for (int j=0; j<qobj2->getSizeVectorStarts(); j++){
                            REQUIRE(qobj2->getVectorStarts()[j]==qstart2[0][j]);
                        }
                    }

    
                }
                SECTION("with linear objectives in the 2nd stage")
                {
                    /** The objective function is:
                     * 150 C0 + 230 C1 + 260 C2 + 79.33333333333333 C3 + 70 C4
                        - 56.66666666666666 C5 - 50 C6 - 12 C7 - 3.333333333333333 C8
                        + 79.33333333333333 C9 + 70 C10 - 56.66666666666666 C11 - 50 C12
                        - 12 C13 - 3.333333333333333 C14 + 79.33333333333333 C15 + 70 C16
                        - 56.66666666666666 C17 - 50 C18 - 12 C19 - 3.333333333333333 C20 + [
                        2 C0 ^2 + 4 C0 * C1 + 10 C1 ^2 ] / 2 
                    */
                    for (int i = 0; i < NS; i++)
                    {
                        loadQuadraticSecondStage(env, i, probability[i], start, index2, value[i], yclbd, ycubd, yctype, profit, NULL, NULL, NULL, NULL, rlbd, rubd);
                    }
                    TssModel *model = getTssModel(env);

                    /** check linear objectives for each scenario in the second stage */
                    const double * elements1;
                    for (int i=0; i<model->getNumScenarios(); i++){
                        elements1=new double [model->getObjScenario(i)->getNumElements()];
                        elements1=model->getObjScenario(i)->getElements();
                        for (int j=0; j<model->getObjScenario(i)->getNumElements(); j++){
                            REQUIRE(elements1[i]==profit[i]);
                        }
                    }
                    
                    /** check quadratic objectives in the first stage */
                    const CoinPackedMatrix *qobj1;
                    qobj1=model->getQuadraticObjCore(0);
                    REQUIRE(qobj1->getMajorDim()==numcol1);
                    for (int i=0; i<qobj1->getNumElements(); i++){
                        REQUIRE(qobj1->getElements()[i]==qvalue[i]);
                        REQUIRE(qobj1->getIndices()[i]==qcolindex[i]);
                    }
                    for (int i=0; i<qobj1->getSizeVectorStarts(); i++){
                        REQUIRE(qobj1->getVectorStarts()[i]==qstart1[i]);
                    }
                    for (int i=0; i<model->getNumScenarios(); i++){
                        REQUIRE(model->getQuadraticObjScenario(i)==NULL);
                    }
                }
            }
            SECTION("with stochastic quadratic objectives in the second stage"){
                    for (int i = 0; i < NS; i++)
                    {
                        /** the objective function is 
                              150 C0 + 230 C1 + 260 C2 + 79.33333333333333 C3 + 70 C4
                                - 56.66666666666666 C5 - 50 C6 - 12 C7 - 3.333333333333333 C8
                                + 79.33333333333333 C9 + 70 C10 - 56.66666666666666 C11 - 50 C12
                                - 12 C13 - 3.333333333333333 C14 + 79.33333333333333 C15 + 70 C16
                                - 56.66666666666666 C17 - 50 C18 - 12 C19 - 3.333333333333333 C20 + [
                                0.6666666666666666 C3 ^2 + 1.333333333333333 C3 * C4
                                + 3.333333333333333 C4 ^2 + 0.6666666666666666 C9 ^2
                                + 0.6666666666666666 C10 ^2 + 1.333333333333333 C11 ^2
                                + 1.333333333333333 C15 ^2 + 0.3333333333333333 C15 * C17
                                + 2.666666666666667 C17 ^2 ] / 2 
                        */
                        loadQuadraticSecondStage(env, i, probability[i], start, index2, value[i], yclbd, ycubd, yctype, profit, qrowindex2[i], qcolindex2[i], qvalue2[i], numq2, rlbd, rubd);
                    }
                    TssModel *model = getTssModel(env);
                    

                    /** check linear objectives for each scenario in the second stage */
                    const double * elements1;
                    for (int i=0; i<model->getNumScenarios(); i++){
                        elements1=new double [model->getObjScenario(i)->getNumElements()];
                        elements1=model->getObjScenario(i)->getElements();
                        for (int j=0; j<model->getObjScenario(i)->getNumElements(); j++){
                            REQUIRE(elements1[i]==profit[i]);
                        }
                    }
                    
                    /** check quadratic objectives in the first stage */
                    const CoinPackedMatrix *qobj1;
                    qobj1=model->getQuadraticObjCore(0);
                    REQUIRE(qobj1->getMajorDim()==numcol1);
                    for (int i=0; i<qobj1->getNumElements(); i++){
                        REQUIRE(qobj1->getElements()[i]==qvalue[i]);
                        REQUIRE(qobj1->getIndices()[i]==qcolindex[i]);
                    }
                    for (int i=0; i<qobj1->getSizeVectorStarts(); i++){
                        REQUIRE(qobj1->getVectorStarts()[i]==qstart1[i]);
                    }

                    /** check quadratic objectives in the second stage */
                    const CoinPackedMatrix *qobj2;
                    
                    for (int i=0; i<model->getNumScenarios(); i++){
                        qobj2=model->getQuadraticObjScenario(i);
                        for (int j=0; j<qobj2->getNumElements(); j++){
                            REQUIRE(qobj2->getElements()[j]==qvalue2[i][j%3]);
                            REQUIRE(qobj2->getIndices()[j]==qcolindex2[i][j%3]);
                        }
        
                        for (int j=0; j<qobj2->getSizeVectorStarts(); j++){
                            REQUIRE(qobj2->getVectorStarts()[j]==qstart2[i][j]);
                        }
                    }
            }
        }
        SECTION("with linear objectives in the 1st stage")
        {
            loadQuadraticFirstStage(env, start1, index1, value1, xclbd, xcubd, xctype, Cost, NULL, NULL, NULL, NULL, rubd1, rubd1);
            TssModel *model = getTssModel(env);
            /** two-stage problem 
              * check number of stages;
              * check number of scenarios;
            */
            REQUIRE(model->getNumStages()==2);
            REQUIRE(model->getNumScenarios()==NS);

            /** number of columns in the first stage problem */
            int numcol1=model->getNumCols(0);
            REQUIRE(numcol1==3);

            /** check linear objectives in each stage */
            /** check first stage core[0] */
            const double *objcore=new double[numcol1];
            objcore=model->getObjCore(0);
            for (int i=0; i<numcol1; i++){
                REQUIRE(objcore[i]==Cost[i]);
            }
            
            SECTION("load second stage")
            {
                SECTION("with quadratic objectives in the 2nd stage")
                {
                    /** The objective function is
                     * 150 C0 + 230 C1 + 260 C2 + 79.33333333333333 C3 + 70 C4
                        - 56.66666666666666 C5 - 50 C6 - 12 C7 - 3.333333333333333 C8
                        + 79.33333333333333 C9 + 70 C10 - 56.66666666666666 C11 - 50 C12
                        - 12 C13 - 3.333333333333333 C14 + 79.33333333333333 C15 + 70 C16
                        - 56.66666666666666 C17 - 50 C18 - 12 C19 - 3.333333333333333 C20 + [
                        0.6666666666666666 C3 ^2 + 1.333333333333333 C3 * C4
                        + 3.333333333333333 C4 ^2 + 0.6666666666666666 C9 ^2
                        + 1.333333333333333 C9 * C10 + 3.333333333333333 C10 ^2
                        + 0.6666666666666666 C15 ^2 + 1.333333333333333 C15 * C16
                        + 3.333333333333333 C16 ^2 ] / 2 
                     */
                    for (int i = 0; i < NS; i++)
                    {
                        loadQuadraticSecondStage(env, i, probability[i], start, index2, value[i], yclbd, ycubd, yctype, profit, qrowindex2[0], qcolindex2[0], qvalue2[0], numq2, rlbd, rubd);
                    }
                    TssModel *model = getTssModel(env);
                    
                    /** check linear objectives for each scenario in the second stage */
                    const double * elements1;
                    for (int i=0; i<model->getNumScenarios(); i++){
                        elements1=new double [model->getObjScenario(i)->getNumElements()];
                        elements1=model->getObjScenario(i)->getElements();
                        for (int j=0; j<model->getObjScenario(i)->getNumElements(); j++){
                            REQUIRE(elements1[i]==profit[i]);
                        }
                    }
                    
                    /** check quadratic objectives in the first stage */
                    
                    REQUIRE(model->getQuadraticObjCore(0)==NULL);

                    /** check quadratic objectives in the second stage */
                    const CoinPackedMatrix *qobj2;
                    
                    for (int i=0; i<model->getNumScenarios(); i++){
                        qobj2=model->getQuadraticObjScenario(i);
                        for (int j=0; j<qobj2->getNumElements(); j++){
                            REQUIRE(qobj2->getElements()[j]==qvalue2[0][j%3]);
                            REQUIRE(qobj2->getIndices()[j]==qcolindex2[0][j%3]);
                        }
        
                        for (int j=0; j<qobj2->getSizeVectorStarts(); j++){
                            REQUIRE(qobj2->getVectorStarts()[j]==qstart2[0][j]);
                        }
                    }
                }
                SECTION("with linear objectives in the 2nd stage")
                {
                    /** This is the original farmers example */
                    for (int i = 0; i < NS; i++)
                    {
                        loadQuadraticSecondStage(env, i, probability[i], start, index2, value[i], yclbd, ycubd, yctype, profit, NULL, NULL, NULL, NULL, rlbd, rubd);
                    }
                    TssModel *model;
                    model = getTssModel(env);
                    /** two-stage problem 
                     * check number of stages;
                     * check number of scenarios;
                    */
                    
                    /** check linear objectives for each scenario in the second stage */
                    const double * elements1;
                    for (int i=0; i<model->getNumScenarios(); i++){
                        elements1=new double [model->getObjScenario(i)->getNumElements()];
                        elements1=model->getObjScenario(i)->getElements();
                        for (int j=0; j<model->getObjScenario(i)->getNumElements(); j++){
                            REQUIRE(elements1[i]==profit[i]);
                        }
                    }
                    
                    /** check quadratic objectives in the first stage */
                    
                    REQUIRE(model->getQuadraticObjCore(0)==NULL);

                    /** check quadratic objectives in the second stage */                    
                    for (int i=0; i<model->getNumScenarios(); i++){
                        REQUIRE(model->getQuadraticObjScenario(i)==NULL);
                    }
                }
            }
        }
    }
}