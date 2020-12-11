#include "catch.hpp"

#include "DspCInterface.h"

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
// const int qrowindex2[3] = {3, 3, 4};
// const int qcolindex2[3] = {3, 4, 4};
// const double qvalue2[3] = {1, 2, 5};
const CoinBigIndex numq2 = 3;
//variable type of second-stage variables
const char yctype[6] = {'C', 'C', 'C', 'C', 'C', 'C'};

//lower and upper limits of constraints: minimum crop requirement
const double rlbd[3] = {200, 240, 0};
const double rubd[3] = {INFINITY, INFINITY, INFINITY};

TEST_CASE("Extensive Formulation for MIQP")
{
    #ifdef DSP_HAS_GRB

    DspApiEnv *env = createEnv();
    /** set number of scenarios */
    setNumberOfScenarios(env, NS);

    /** set problem dimension */
    setDimensions(env, 3, 1, 6, 3);

    SECTION("load first stage")
    {
        SECTION("with quadratic objectives in the 1st stage")
        {
            loadQuadraticFirstStage(env, start1, index1, value1, xclbd, xcubd, xctype, Cost, qrowindex, qcolindex, qvalue, numq, rubd1, rubd1);
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
                    /** using Gurobi to solve extensive formulation */
                    setIntParam(env, "DE/SOLVER", 4);
                    solveDe(env);
                    double prim;
                    prim = getPrimalBound(env);
                    double p = round(prim * 1000000) / 1000000;
                    REQUIRE(p == -42666.988254);
                    //free memory
                    freeSolver(env);
                    freeModel(env);
                    freeEnv(env);
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
                    setIntParam(env, "DE/SOLVER", 4);
                    solveDe(env);
                    double prim;
                    prim = getPrimalBound(env);
                    double p = round(prim * 1000000) / 1000000;
                    REQUIRE(p == -68826.562500);
                    //free memory
                    freeSolver(env);
                    freeModel(env);
                    freeEnv(env);
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
                    /** using Gurobi to solve extensive formulation */
                    setIntParam(env, "DE/SOLVER", 4);
                    solveDe(env);
                    double prim;
                    prim = getPrimalBound(env);
                    double p = round(prim * 1000000) / 1000000;
                    REQUIRE(p == -55948.852793);
                    //free memory
                    freeSolver(env);
                    freeModel(env);
                    freeEnv(env);
            }
        }
        SECTION("with linear objectives in the 1st stage")
        {
            loadQuadraticFirstStage(env, start1, index1, value1, xclbd, xcubd, xctype, Cost, NULL, NULL, NULL, NULL, rubd1, rubd1);
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
                    setIntParam(env, "DE/SOLVER", 4);
                    solveDe(env);
                    double prim;
                    prim = getPrimalBound(env);
                    double p = round(prim * 1000000) / 1000000;
                    REQUIRE(p == -108251.276042);
                    //free memory
                    freeSolver(env);
                    freeModel(env);
                    freeEnv(env);
                }
                SECTION("with linear objectives")
                {
                    /** This is the original farmers example */
                    for (int i = 0; i < NS; i++)
                    {
                        loadQuadraticSecondStage(env, i, probability[i], start, index2, value[i], yclbd, ycubd, yctype, profit, NULL, NULL, NULL, NULL, rlbd, rubd);
                    }
                    setIntParam(env, "DE/SOLVER", 4);
                    solveDe(env);
                    double prim;
                    prim = getPrimalBound(env);
                    double p = round(prim * 1000000) / 1000000;
                    REQUIRE(p == -108390.000000);
                    //free memory
                    freeSolver(env);
                    freeModel(env);
                    freeEnv(env);
                }
            }
        }
    }
    #endif
}