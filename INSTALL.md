Installation
============
If you have all the prerequisite packages installed on your system (see the "BUILD ESSENTIALS" and "PREREQUISITES" sections below), then you need to go to the root directory of DSP and type
```cmake
cmake .
```
to configure OOQP. If you wish to install the package in a more permanent location, you may then type
```
make install
```
External packages (MA27, OOQP, SCIP Optimization Suite, Smi) used in DSP are built automatically. A shared object is installed in ./lib directory. Once the installation has been successfully done, you need to set an environment variable DSP_INC and (DY)LD_LIBRARY_PATH.

For Linux,
```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<DSP source path>/lib
```
For Mac,
```
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:<DSP source path>/lib
```

BUILD ESSENTIALS
----------------
The following packages are essential to build DSP and the required external software packages. If you do not have any of the following packages, and if you do not know how to get it, please see the "HOW TO GET BUILD ESSENTIALS" section below.
* CMake
* GNU Autoconf
* GNU Automake
* GNU Make
* MPICH
* BLAS
* LAPACK
* SVN

PREREQUISITES
-------------
The following packages need to be installed on your system and located on ./extra directory before DSP may be built.

* MA27 (./extra/ma27-1.0.0) -- MA27 is a library for solving sparse symmetric indefinite linear systems. To build OOQP solver, you must have this installed. This is part of HSL (formerly the Harwell Subroutine Library), a collection of ISO Fortran codes for large scale scientific computation. You can download the source files and obtain the free personal license from http://hsl.rl.ac.uk/archive/hslarchive.html. Move the downloaded directory to ./extras/ma27-1.0.0/

* OOQP (./extra/OOQP-0.99.25) -- OOQP is an object-oriented C++ package, based on a primal-dual interior-point method, for solving convex quadratic programming problems. You can get the source code per request to the developer: http://pages.cs.wisc.edu/~swright/ooqp/

* SCIP Optimization Suite (./extra/scipoptsuite-3.1.1) -- SCIP is non-commercial solvers for mixed integer programming (MIP) and mixed-integer nonlinear programming (MINLP). It is freely available from the website: http://scip.zib.de/download.php?fname=scipoptsuite-3.1.1.tgz

JULIA
-----
Julia is a high level dynamic programming language for technical computing, with syntax that is familiar to users of other technical computing environments such as MATLAB and Python. DSP provides an interface to an algebraic modeling package StochJuMP based on Julia. The interface package can be installed by the Julia command
```julia
Pkg.clone("https://github.com/kibaekkim/DSP.jl.git");
```
and the examples are located at
```
./examples/julia
```
Julia can be downloaded from http://julialang.org/downloads/
