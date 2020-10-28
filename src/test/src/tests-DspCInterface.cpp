// tests-BlkModel.cpp
#include "catch.hpp"

#include "DspCInterface.h"

TEST_CASE("C interface functions") {
    DspApiEnv* env = NULL;

    SECTION("createEnv function") {
        env = createEnv();
        REQUIRE(env != NULL);
        REQUIRE(env->model_ == NULL);
        REQUIRE(env->solver_ == NULL);
    }

    int NS = 3;                        // number of scenarios
    double probability[3] = {(double)1/3, (double)1/3, (double)1/3}; // probability

    //first-stage constraint matrix
    CoinBigIndex start1[2]={0,3};
    const int index1[3]={0, 1, 2};
    const double value1[3]={1, 1, 1};

    //quadratic objective information
    const int qrowindex[3]={0, 0, 1};
    const int qcolindex[3]={0, 1, 1};
    const double qvalue[3]={1, 2, 5};
    const CoinBigIndex numq=3;

    //first-stage variable lower bound, upper bound, and type
    const double xclbd[3]={0, 0, 0};
    const double xcubd[3]={500, 500, 500};
    const char xctype[3]={'C', 'C', 'C'};

    //first-stage objective 
    const double Cost[3] = {150, 230, 260};   // cost of planting crops

    //first-stage constraint upper and lower limit
    const double rubd1[1]   = {500};             // budget capacity
    const double rlbd1[1]   = {0};             // budget capacity

    

    /** second-stage constraint matrix
     *  the index starts from the first-stage variables,
     *  and the maximum is less than 
     *  # first-stage variables + # of single scenario second-stage variables
     **/
    //second-stage objective coefficient: only containes second stage variables
    //selling price 238 210
    //purchase price -170, -150, -36, -10
    const double profit[6]= {238, 210, -170, -150, -36, -10};

    //second-stage constraint matrix
    const CoinBigIndex start[4]={0, 3, 6, 9};
    const int index[9]={0, 3, 5, 1, 4, 6, 2, 7, 8};
    const double value[3][9]={
        {3.0, 1, -1, 3.6, 1, -1, 24.0, -1, -1},
        {2.5, 1, -1, 3.0, 1, -1, 20.0, -1, -1},
        {2.0, 1, -1, 2.4, 1, -1, 16.0, -1, -1}
    };

    //lower and upper limits of second-stage variables
    const double yclbd[6]={0, 0, 0, 0, 0, 0};
    const double ycubd[6]={INFINITY, INFINITY, INFINITY, INFINITY, 6000, INFINITY};

    //second stage quadratic objective
    const int qrowindex2[3]={3, 3, 4};
    const int qcolindex2[3]={3, 4, 4};
    const double qvalue2[3]={1, 2, 5};
    const CoinBigIndex numq2=3;
    //variable type of second-stage variables
    const char yctype[6]={'C', 'C', 'C', 'C', 'C', 'C'};

    //lower and upper limits of constraints: minimum crop requirement
    const double rlbd[3]={200, 240, 0};
    const double rubd[3]={INFINITY, INFINITY, INFINITY};

    /** set number of scenarios */
    //setNumberOfScenarios(env, NS);

    SECTION("load quadratic first stage problem"){
        //env = createEnv();
        //setNumberOfScenarios(env, NS);
        //setDimensions(env, 3, 1, 6, 3);
        //loadQuadraticFirstStage(env, start1, index1, value1, xclbd, xcubd, xctype, Cost, qrowindex, qcolindex, qvalue, numq, rubd1, rubd1);
        //for (int i=0; i<NS; i++){
        //    loadQuadraticSecondStage(env, i, probability[i], start, index, value[i], yclbd, ycubd, yctype, profit,qrowindex2, qcolindex2, qvalue2, numq2, rlbd, rubd);
        //}
        //setIntParam(env, "DE/SOLVER", 4);
        //solveDe(env);
        //double prim;
        //prim=getPrimalBound(env);
        //double p=round(prim*1000000)/1000000;
        //REQUIRE(p == -64818.588589);
    }

    SECTION("freeEnv function") {
        freeEnv(env);
        REQUIRE(env == NULL);
    }
}