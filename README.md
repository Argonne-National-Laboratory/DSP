DSP
===
DSP is an open-source optimization solver that implements decomposition methods for stochastic mixed-integer programming (SMIP) problems. These are structured optimization problems considering uncertain scenario realization s with probability p_s in the following form:

        minimize   c^T x + \sum_{s=1}^S p_s q_s^T y_s
        subject to   A x                              = b
                   T_s x +                    W_s y_s = h_s for s = 1, .., S
                   some x, y_s are integers

where x and y_s are decision variable vectors with dimensions n_1 and n_2, respectively, A, T_s and W_s are matrices of dimensions m_1 by n_1, m_2 by n_1 and m_2 by n_2, respectively, and c, q_s, b, and h_s are vectors of appropriate dimensions.

The problem structures allow decomposition approaches for solving the problem in parallel computing system. DSP provides parallel implementations for decomposition methods (Benders decomposition and dual decomposition). The methods can be run on cluster as well as desktop computers.

Prerequisites
-------------
We recommend to install and run DSP on a linux-type cluster with MPI library. If a cluster is not available for you, then the following packages are required to build DSP. Please have the packages available on your machine.
* MPICH -- This is **required** to build and run DSP. It is required to install a version of MPICH, which is available from http://www.mpich.org/downloads/. On Linux, you can also do
```bash
    sudo apt-get install libmpich-dev
```
* CMake -- This package is **required** to build DSP and all the other external software packages. It is available from http://www.cmake.org/download/. On Linux, you can also do
```bash
    sudo apt-get install cmake
```
* BLAS/LAPACK -- These linear algebra libraries are **required** and may already be available on your machine. If not, the source codes are available from http://www.netlib.org/blas/blas.tgz and http://www.netlib.org/lapack/lapack-3.5.0.tgz. Or on linux you can install them by typing
```bash
    sudo apt-get install libblas-dev
    sudo apt-get install liblapack-dev
```
* svn -- A subversion software package is **required** to download an external package Coin-SMI used in DSP. You can install it by
```bash
    sudo apt-get install subversion
```
* GNU Make/Autoconf/Automake -- These are **required** and likely available on any machine. If not, you can get them as follows.
```bash
    sudo apt-get install build-essential
    sudo apt-get install autoconf
    sudo apt-get install automake
```

The following packages are also **required** to build and run DSP and need to be located on ./extra directory before DSP may be built. **DSP will automatically configure and build the external packages once they are located in the right place.**

* MA27 (./extra/ma27-1.0.0) -- MA27 is a library for solving sparse symmetric indefinite linear systems. To build OOQP solver, you must have this installed. This is part of HSL (formerly the Harwell Subroutine Library), a collection of ISO Fortran codes for large scale scientific computation. You can download the source files and obtain a free personal license from http://www.hsl.rl.ac.uk/download/MA27/1.0.0/a/. Move the downloaded directory to ./extras/ma27-1.0.0/

* SCIP Optimization Suite (./extra/scipoptsuite-3.1.1) -- SCIP is non-commercial solvers for mixed integer programming (MIP) and mixed-integer nonlinear programming (MINLP). It is freely available from the website: http://scip.zib.de/download.php?fname=scipoptsuite-3.1.1.tgz

Installation
------------
If you have all the prerequisite packages installed on your system, then you need to go to the root directory of DSP and type
```bash
cmake .
```
to configure OOQP. If you wish to install the package in a more permanent location, you may then type
```bash
make install
```
External packages (MA27, OOQP, SCIP Optimization Suite, Smi) used in DSP are built automatically. A shared object is installed in ./lib directory. Once the installation has been successfully done, you need to set environment variable (DY)LD_LIBRARY_PATH.
For Linux,
```bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<DSP_SRC_PATH>/lib
```
For Mac,
```bash
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:<DSP_SRC_PATH>/lib
```

Julia interface (Required)
--------------------------
DSP uses Julia as a modeling interface for the computatioanl experiments. Julia is a high level dynamic programming language for technical computing, with syntax that is familiar to users of other technical computing environments such as MATLAB and Python. Julia can be downloaded from http://julialang.org/downloads/ If Julia is successfully installed on your machine, then you can start the Julia command-line tool by typing
```bash
julia
```
Now we need to install four Julia packages required to run DSP in the Julia environment. The packages should be installed in the Julia command-line tool.
* MPI.jl is an MPI interface package, which can be installed by the following Julia command. MPICH is required for this package.
```julia
Pkg.add("MPI");
```
* JuMP.jl is a algebraic modeling package in Julia for mathematical programming, which can be installed by the following Julia command:
```julia
Pkg.add("JuMP");
```
* StochJuMP.jl is a algebraic modeling package in Julia for stochastic programming, which can be installed by the following Julia command:
```julia
Pkg.clone("https://github.com/kibaekkim/StochJuMP.jl.git");
```
* The DSPsolver.jl package provides an interface to StochJuMP. DSPsolver.jl can be installed by the Julia command
```julia
Pkg.clone("https://github.com/kibaekkim/DSPsolver.jl.git");
```

Reproducing the Computational Results
-------------------------------------
We provide scripts to reproduce the computational results reported in the manuscript. The scripts are available in subdirectory ./experiments/ that contains the bash-shell scripts for running DSP and the MATLAB scripts for reproducing tables and figures.
* Table2.pbs.sh -- This is a bash script that runs the computational experiments for Table 2 in the manuscript on cluster with the PBS scheduler. If a cluster is not available, then Table2.sh can be used on local machine. The results will be stored in subdirectory ./output in CSV format.
 * Figures2and3.m -- This MATLAB script will generate Figure2.eps and Figure3.eps from the CSV files generated by Table2.pbs.sh (or Table2.sh).
 * Table2.m -- This MATLAB script will generate Table2.tex from the CSV files generated by Table2.pbs.sh (or Table2.sh).
 * Figure4.m -- This MATLAB script will generate Figure4.eps from the CSV files generated by Table2.pbs.sh (or Table2.sh).
* Figure5.pbs.sh
* Figure6.pbs.sh
* Table4.pbs.sh
* Figure8.pbs.sh
* Table5.pbs.sh

Contents of the distribution
----------------------------
The DSP distribution contains the following top-level subdirectories:
* experiments/ -- All the scripts to reproduce computational results including tables and figures in the manuscript.
* extra/ -- External packages required and used for DSP. These include MA27, OOQP, SCIP and SMI.
* farmer/ -- This contains an example for using DSP with Julia.
* lib/ -- This contains DSP shared object precompiled on a 64-bit linux machine.
* src/ -- Source files for the DSP distribution.

Credits
-------
DSP has been developed and is maintained by:
* Kibaek Kim, Mathematics and Computer Science Division, Argonne National Laboratory.
* Victor M. Zavala, Mathematics and Computer Science Division, Argonne National Laboratory.

Acknowledgements
----------------
We thank E. Michael Gertz and Stephen Wright for providing OOQP software package.
