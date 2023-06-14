#include <stdlib.h>
#include <stdio.h>
#include <dlfcn.h>

/*
Farmer example: This is based on a well-known stochastic programming problem.

minimize 1/3 * (150 x11 + 230 x21 + 260 x31 + 238 y11 + 210 y21 - 170 w11 - 150 w21 - 36 w31 - 10 w41)
       + 1/3 * (150 x12 + 230 x22 + 260 x32 + 238 y12 + 210 y22 - 170 w12 - 150 w22 - 36 w32 - 10 w42)
       + 1/3 * (150 x13 + 230 x23 + 260 x33 + 238 y13 + 210 y23 - 170 w13 - 150 w23 - 36 w33 - 10 w43)
subject to
  x11             - x12                               == 0
        x21             - x22                         == 0
              x31             - x32                   == 0
                    x12             - x13             == 0
                          x22             - x23       == 0
                                x32             - x33 == 0
  x11 + x21 + x31                                     <= 500
  
Each block has the following constraint matrix and the common row bound vectors:

    A1 = [3.0         1.0     -1.0               ]
         [    3.6         1.0      -1.0          ]
         [         24                   -1.0 -1.0]
         [                               1.0     ]

    A2 = [2.5         1.0     -1.0               ]
         [    3.0         1.0      -1.0          ]
         [         20                   -1.0 -1.0]
         [                               1.0     ]

    A3 = [2.0         1.0     -1.0               ]
         [    2.4         1.0      -1.0          ]
         [         16                   -1.0 -1.0]
         [                               1.0     ]

    bl = [200 240 0.0 -Inf]
    bu = [Inf Inf Inf 6000]
*/

int main(int argc, char **argv)
{
    void *handle;
    char *error;

    // Load the library
    handle = dlopen("libDsp.dylib", RTLD_LAZY);
    if (!handle) {
        fputs (dlerror(), stderr);
        exit(1);
    }

    void* env = NULL;
    env = ((void *(*)())dlsym(handle, "createEnv"))();
    if (env == NULL) {
        printf("Failed to create DSP environment.\n");
    } else {
        printf("Successfully created DSP environment.\n");
    }

    // Problem Data
    int nvars = 27;
    int ncons0 = 7;
    int ncons1 = 4;
    double obj[] = {150/3, 230/3, 260/3, 238/3, 210/3, -170/3, -150/3, -36/3, -10/3};
    double rlbd0[] = {0, 0, 0, 0, 0, 0, -1e+21};
    double rubd0[] = {0, 0, 0, 0, 0, 0, 500};
    double rlbd1[] = {200, 240, 0, -1e+21};
    double rubd1[] = {1e+21, 1e+21, 1e+21, 6000};
    int start0[] = {0, 2, 4, 6, 8, 10, 12, 15};
    int index0[] = {0, 9, 1, 10, 2, 11, 9, 18, 10, 19, 11, 20, 0, 1, 2};
    double value0[] = {1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, 1, 1};
    int start1[] = {0, 3, 6, 9, 10};
    int index1[] = {0, 3, 5, 1, 4, 6, 2, 7, 8, 7};
    double value1[] = {3, 1, -1, 3.6, 1, -1, 24, -1, -1, 1};
    double value2[] = {2.5, 1, -1, 3, 1, -1, 20, -1, -1, 1};
    double value3[] = {2, 1, -1, 2.4, 1, -1, 16, -1, -1, 1};

    double *obj0 = new double [nvars];
    double *clbd0 = new double [nvars];
    double *cubd0 = new double [nvars];
    char *ctype0 = new char [nvars];
    for (int j = 0; j < nvars; ++j) {
        obj0[j] = obj[j % 9];
        clbd0[j] = 0.0;
        cubd0[j] = 1e+21;
        ctype0[j] = j % 9 < 3 ? 'I' : 'C';
    }

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
        nvars,
        ncons0,
        start0[ncons0],
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

    int *index_ = new int [start1[ncons1]];
    double *value_ = new double [start1[ncons1]];

    // create each block
    for (int block = 0; block < 3; block++) {
        for (int j = 0; j < start1[ncons1]; ++j) {
            index_[j] = index1[j] + block * 9;
            if (block == 0) value_[j] = value1[j];
            else if (block == 1) value_[j] = value2[j];
            else if (block == 2) value_[j] = value3[j];
        }
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
            block + 1,
            nvars,
            ncons1,
            start1[ncons1],
            start1,
            index_,
            value_,
            clbd0,
            cubd0,
            ctype0,
            obj0,
            rlbd1,
            rubd1
        );
    }

    // Important to call this!
    ((void (*)(void*))dlsym(handle, "updateBlocks"))(env);

    // Solve the deterministic equivalent problem
    // ((void (*)(void*))dlsym(handle, "solveDe"))(env);

    // FIXME: This causes segfault, probably due to incorrectly addressing unbounded subproblem solution.
    ((void (*)(void*))dlsym(handle, "solveDw"))(env);

    double primal_bound = ((double (*)(void*))dlsym(handle, "getPrimalBound"))(env);
    printf("Primal bound = %e\n", primal_bound);

    double *solution = new double [nvars];
    ((void (*)(void*, int, double*))dlsym(handle, "getPrimalSolution"))(env, nvars, solution);
    for (int j = 0; j < nvars; ++j) {
        printf("var%2d = %f\n", j, solution[j]);
    }
    // ((void (*)(void*, const char*))dlsym(handle, "writeMps"))(env, "farmer");
    
    if (env != NULL) {
        ((void (*)(void*&))dlsym(handle, "freeEnv"))(env);
        env = NULL;
        printf("Successfully freed DSP environment.\n");
    }
    dlclose(handle);

    delete [] obj0;
    delete [] clbd0;
    delete [] cubd0;
    delete [] ctype0;
    delete [] index_;
    delete [] value_;
    delete [] solution;

    return 0;
}
