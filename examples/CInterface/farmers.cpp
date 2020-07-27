#include <iostream>
#include "../src/DspCInterface.h"

/** This is the example of Farmers problem in Section 1.1, Intorduction to Stochatic Programming 
 *  Refer to DspCInterface.hpp for details of functions.
*/

int main(){
    /** create environment variable */
    DspApiEnv* env = createEnv();

    int NS = 3;                        // number of scenarios
    double probability[3] = {(double)1/3, (double)1/3, (double)1/3}; // probability

    /** first-stage constraint matrix
     *  x1+x2+x3<=500 */
    CoinBigIndex start1[1]={0};
    const int index1[3]={0, 1, 2};
    const double value1[3]={1, 1, 1};

    /** first-stage variable lower bound, upper bound, and type */
    const double xclbd[3]={0, 0, 0};
    const double xcubd[3]={INFINITY, INFINITY, INFINITY};
    const char xctype[3]={'C', 'C', 'C'};

    /** first-stage objective */
    const double Cost[3] = {150, 230, 260};   // cost of planting crops

    /** first-stage constraint upper and lower limit */
    const double rubd1[1]   = {500};             // budget capacity
    const double rlbd1[1]   = {0};             // budget capacity

    /** second-stage constraint matrix
     *  the index starts from the first-stage variables,
     *  and the maximum is less than 
     *  # first-stage variables + # of second-stage variables in a single scenario
     **/

    /** second-stage objective coefficient: only containes second stage variables
     *  selling price 238 210
     *  purchase price -170, -150, -36, -10*/
    const double profit[6]= {238, 210, -170, -150, -36, -10};

    /** second-stage constraint matrix*/
    const CoinBigIndex start[4]={0, 3, 6, 9};
    const int index[9]={0, 3, 5, 1, 4, 6, 2, 7, 8};
    const double value[3][9]={
        {3.0, 1, -1, 3.6, 1, -1, 24.0, -1, -1},
        {2.5, 1, -1, 3.0, 1, -1, 20.0, -1, -1},
        {2.0, 1, -1, 2.4, 1, -1, 16.0, -1, -1}
    };

    /** lower and upper limits of second-stage variables */
    const double yclbd[6]={0, 0, 0, 0, 0, 0};
    const double ycubd[6]={INFINITY, INFINITY, INFINITY, INFINITY, 6000, INFINITY};

    /** variable type of second-stage variables */
    const char yctype[6]={'C', 'C', 'C', 'C', 'C', 'C'};

    /** lower and upper limits of constraints: minimum crop requirement */
    const double rlbd2[3]={200, 240, 0};
    const double rubd2[3]={INFINITY, INFINITY, INFINITY};

    /** set number of scenarios */
    setNumberOfScenarios(env, NS);

    /** set problem dimension */
    setDimensions(env, 3, 1, 6, 3);
    
    /** load first-stage problem */
    loadFirstStage(env, start1, index1, value1, xclbd, xcubd, xctype, Cost, rubd1, rubd1);
    
    /** load second-stage problem corresponding to each sceanrios */
    for (int i=0; i<NS; i++){
        loadSecondStage(env, i, probability[i], start, index, value[i], yclbd, ycubd, yctype, profit, rlbd2, rubd2);
    }

    /** write out MPS files */
    writeMps(env, "farmer");

    /** set solver for deterministic equivalent model, "0" indicates using CPLEX */
    setIntParam(env, "DE/SOLVER", 0);

    /** solve deterministic equivalent model */
    solveDe(env);

    /** print primal solution */
    double prim;
    prim=getPrimalBound(env);
    printf("optimal value = %f", prim);

    /** free memory */
    freeSolver(env);
    freeModel(env);
    freeEnv(env);

}
