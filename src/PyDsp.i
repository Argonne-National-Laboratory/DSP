%module PyDsp
%{
#include "StoCInterface.h"
%}

%include "carrays.i"
%array_class(double, doubleArray);

%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);

StoApiEnv * createEnv();
void freeEnv(StoApiEnv * env);

void readSmps(StoApiEnv * env, const char * smps);

void readParamFile(StoApiEnv * env, const char * param_file);

void solveDe(StoApiEnv * env);
void solveBd(StoApiEnv * env, int nauxvars);
void solveBdMpi(StoApiEnv * env, int nauxvars, MPI_Comm comm);
void solveDd(StoApiEnv * env, MPI_Comm comm);

int getSolutionStatus(StoApiEnv * env);
int getTotalNumCols(StoApiEnv * env);
double getPrimalBound(StoApiEnv * env);
double getDualBound(StoApiEnv * env);
void getSolution(StoApiEnv * env, int num, double * solution);
