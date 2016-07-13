%module PyDsp
%{
#include "DspCInterface.h"
%}

%include "carrays.i"
%array_class(double, doubleArray);

%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);

DspApiEnv * createEnv();
void freeEnv(DspApiEnv * env);

void readSmps(DspApiEnv * env, const char * smps);

void readParamFile(DspApiEnv * env, const char * param_file);

void solveDe(DspApiEnv * env);
void solveDd(DspApiEnv * env);
void solveBd(DspApiEnv * env, int nauxvars);
void solveDdMpi(DspApiEnv * env, MPI_Comm comm);
void solveBdMpi(DspApiEnv * env, int nauxvars, MPI_Comm comm);

int getStatus(DspApiEnv * env);
int getTotalNumCols(DspApiEnv * env);
double getPrimalBound(DspApiEnv * env);
double getDualBound(DspApiEnv * env);
void getSolution(DspApiEnv * env, int num, double * solution);
