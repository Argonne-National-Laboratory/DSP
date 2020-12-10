# Installation

`DSP` is an open source and available in a Github repository.
We describe

* how to correctly clone the repository with submodules
* how to install DSP and Julia interface
* the prerequisites for installation

## Download

You can clone this repository in your preferred directory by typing

=== "Single line"

    ```
    git clone --recursive https://github.com/Argonne-National-Laboratory/DSP.git
    ```

=== "Multiple lines"

    ```
    git clone https://github.com/Argonne-National-Laboratory/DSP.git
    cd DSP
    git submodule update --init --recursive
    ```

## Build and Setting

Please set paths in `UserConfig.cmake` as follows.

```
set(MA27LIB_DIR     "/path/to/libma27/dir")
set(CPLEX_LIB_DIR   "/path/to/libcplex/dir")
set(CPLEX_INC_DIR   "/path/to/ilcplex")
set(GUROBI_LIB_DIR  "/path/to/libgrb/dir")
set(GUROBI_INC_DIR  "/path/to/gurobi/include")
set(SCIPOPT_INC_DIR "/path/to/scip/include")
set(SCIPOPT_LIB_DIR "/path/to/scip/lib")
```

!!! info
    Not every path is required to set in `UserConfig.cmake`. User needs to set whatever they need.

    For example, `MA27LIB_DIR` is required only if user wants to use `OOQP` solver for the Lagrangian master problem in dual decomposition.
    And, most users are likely to need only one of the solvers: `CPLEX`, `GUROBI`, or `SCIP`.

Assuming that you are at the root directory of DSP, type

```
mkdir build
cd build
cmake ..
```

If everything goes well, you should be able to find these two files: `./bin/runDsp` and `./lib/libDsp.*`.
Optionally, you can test the build with `ctest`.

---

## Requirement

We recommend to install and run DSP on **Mac** or **Linux** machine optionally with an appropriate **MPI** library. 
The software packages necessary to build the source code of DSP are listed below.

### External Solvers

At least one of the following external MIP solvers needs to be given in `UserConfig.cmake`.

- CPLEX
- Gurobi
- SCIP: We strongly recomment to use SCIP version 7 or higher.

OOQP can be optionally used by setting the path to MA27 library.a

### Libraries

The following libraries are necessary but may already be available on most computers.
If not, you can download them by using `apt-get`, `brew`, or any other way of your preference.

- GNU/LLVM C++ compiler: the other compilers (i.e., intel) were not tested.
- [CMake](http://www.cmake.org/download/): This package is to build DSP and all the other external software packages.
- BLAS/LAPACK: These linear algebra libraries may already be available on your machine.
- gfortran: GNU Fortran compiler and library.
- bzip2/zlib: These packages are required to build one of the external software packages used in DSP.

### Optional Library

- MPI: This library is optional to build and run DSP, but required to run DSP in parallel.