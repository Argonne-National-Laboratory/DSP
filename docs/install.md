# Installation

This describes the following:

* how to correctly clone the repository with submodules
* the prerequisites for installation
* how to install DSP and Julia interface

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

Please set `UserConfig.cmake` as follows.
Assuming that you are at the root directory of DSP, type

```
mkdir build
cd build
cmake ..
```

A shared object is installed in ``./lib`` directory. 
Once the installation has been successfully done, you need to set environment variable ``(DY)LD_LIBRARY_PATH``.

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