Prerequisites
-------------

We recommend to install and run DSP on a **Linux** machine with an appropriate
**MPI** library. The software packages necessary to build the source code of
DSP are listed below.

Required Libraries
^^^^^^^^^^^^^^^^^^

* CMake -- This package is to build DSP and all the other external software packages. It is available from http://www.cmake.org/download/.
* BLAS/LAPACK -- These linear algebra libraries may already be available on your machine. If not available, the source codes are available from http://www.netlib.org/blas/blas.tgz and http://www.netlib.org/lapack/lapack-3.5.0.tgz.
* GNU/LLVM C++ compiler -- the other compilers (i.e., intel) were not tested.
* gfortran -- GNU Fortran compiler and library.
* bzip2/zlib -- These packages are required to build one of the external software packages used in DSP.

Optional Library
^^^^^^^^^^^^^^^^

* MPI - This library is optional to build and run DSP, but required to run DSP in parallel.

External Optimization Solvers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* MA27
* Cplex
* SCIP