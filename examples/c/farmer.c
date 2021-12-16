#include <stdlib.h>
#include <stdio.h>
#include <dlfcn.h>

int main(int argc, char **argv)
{
    void *handle;
    char *error;

    handle = dlopen ("libDsp.dylib", RTLD_LAZY);
    if (!handle) {
        fputs (dlerror(), stderr);
        exit(1);
    }

    void* env = NULL;
    env = ((void *(*)())dlsym(handle, "createEnv"))();
    if (env == NULL) {
        printf("Failed to create DSP environment.\n");
    } else {
        printf("Successful to create DSP environment.\n");
    }

    // Problem Data
    int NS = 3;
    double probability[] = {1/3, 1/3, 1/3};
    int CROPS = 2;
    int PURCH = 3;
    int SELL = 4;
    double Cost[] = {150, 230, 260};
    double Budget = 500;
    double Purchase[] = {238, 210};
    double Sell[] = {170, 150, 36, 10};
    double Yield[] = {3.0, 3.6, 24.0, 2.5, 3.0, 20.0, 2.0, 2.4, 16.0};
    double Minreq[] = {200, 240, 0};

    /*
    Formulation

    vars: x[1,1], x[2,1], ..., y[1,1], y[2,1], ..., w[1,1], w[2,1], ...

    A0 = [
        1 0 -1  0 ...
        0 1  0 -1 ...
        0 0  1  0 -1  0 ...
        0 0  0  1  0 -1 ...
    ]
    */

    // number of variables
    int nx = CROPS * NS;
    int ny = PURCH * NS;
    int nw = SELL * NS;
    int nvars0 = nx + ny + nw;

    // number of constraints
    int ncons0 = CROPS * (NS - 1);
    int nnz0 = ncons0 * 2;
    int start0[] = {0, 2, 4, 6, 8};
    int index0[] = {0, 2, 1, 3, 2, 4, 3, 5};
    double value0[] = {1, -1, 1, -1, 1, -1, 1, -1};
    double *clbd0 = new double [nvars0];
    double *cubd0 = new double [nvars0];
    char *ctype0 = new char [nvars0];
    double *obj0 = new double [nvars0];
    for (int j = 0; j < nvars0; ++j) {
        clbd0[j] = 0.0;
        cubd0[j] = 1.0e+20;
        ctype0[j] = j < nx ? 'I' : 'C';
    }
    for (int s = 0; s < NS; ++s) {
        for (int i = 0; i < CROPS; ++i)
            obj0[s*CROPS+i] = probability[s] * Cost[i];
        for (int i = 0; i < PURCH; ++i)
            obj0[nx+s*PURCH+i] = probability[s] * Purchase[i];
        for (int i = 0; i < SELL; ++i)
            obj0[nx+nw+s*SELL+i] = probability[s] * Sell[i];
    }
    double rlbd0[] = {0.0, 0.0, 0.0, 0.0};
    double rubd0[] = {0.0, 0.0, 0.0, 0.0};

    ((void (*)(
        void *,
        int,
        int,
        int,
        int,
        int*,
        int*,
        double*,
        double*,
        double*,
        char*,
        double*,
        double*,
        double*
    ))dlsym(handle, "loadBlockProblem"))(
        env,
        0,
        nvars0,
        ncons0,
        nnz0,
        start0,
        index0,
        value0,
        clbd0,
        cubd0,
        ctype0,
        obj0,
        rlbd0,
        rubd0
    );

    /*
    First block subproblem

    A1 = [
        ... 1 1 1 0 ...
        ... x 0 0 1 0 -1  0  0  0 ...
        ... 0 x 0 0 1  0 -1  0  0 ...
        ... 0 0 x 0 0  0  0 -1 -1 ...
        ... 0 0 0 0 0  0  0  1  0 ...
    ]
    */

    int ncons1 = PURCH + 3;
    int nnz1 = CROPS + 3 * PURCH + 3 + 1;
    int start1[] = {};
    int index1[] = {};

    ((void (*)(
        void *,
        int,
        int,
        int,
        int,
        int*,
        int*,
        double*,
        double*,
        double*,
        char*,
        double*,
        double*,
        double*
    ))dlsym(handle, "loadBlockProblem"))(
        env,
        1,
        nvars0,
        ncons1,
        nnz1,
        start1,
        index1,
        value1,
        clbd0,
        cubd0,
        ctype0,
        obj0,
        rlbd1,
        rubd1
    );
    
    if (env != NULL) {
        ((void (*)(void*&))dlsym(handle, "freeEnv"))(env);
        env = NULL;
    }
    dlclose(handle);

    return 0;
}
