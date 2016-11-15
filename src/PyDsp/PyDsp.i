%module PyDsp
%{
#include "DspCInterface.cpp"
%}

%include "carrays.i"
%array_class(double, doubleArray);

DspApiEnv * createEnv();
void freeEnv(DspApiEnv * env);

void readSmps(DspApiEnv * env, const char * smps);

void readParamFile(DspApiEnv * env, const char * param_file);

void solveDe(DspApiEnv * env);
void solveDd(DspApiEnv * env);
void solveBd(DspApiEnv * env);

int getStatus(DspApiEnv * env);
int getTotalNumCols(DspApiEnv * env);
double getPrimalBound(DspApiEnv * env);
double getDualBound(DspApiEnv * env);
void getPrimalSolution(DspApiEnv * env, int num, double * solution);

void writeMps(DspApiEnv * env, const char * name);
