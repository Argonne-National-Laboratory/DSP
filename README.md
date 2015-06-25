# DSP
DSP is an open-source package that implements decomposition methods for **stochastic mixed-integer programming (SMIP)** problems. These are structured optimization problems considering uncertain scenario realizations s with probabilities p_s in the following form:

        minimize   c^T x + \sum_{s=1}^S p_s q_s^T y_s
        subject to   A x                              = b
                   T_s x +                    W_s y_s = h_s for s = 1, .., S
                   some x, y_s are integers

where x and y_s are decision variable vectors with dimensions n_1 and n_2, respectively, A, T_s and W_s are matrices of dimensions m_1 by n_1, m_2 by n_1 and m_2 by n_2, respectively, and c, q_s, b, and h_s are vectors of appropriate dimensions.

DSP provides **parallel** implementations for decomposition methods:
* Dual decomposition
* Benders decomposition

The methods can be run on computing clusters and multi-core processors.

## Download

You can clone this repository in your preferred directory by typing
```bash
git clone https://github.com/kibaekkim/DSP.git
```
or download the zip file from [this link](https://github.com/kibaekkim/DSP/archive/MPC.zip).

## Prerequisites

We recommend to install and run DSP on a **Linux** machine with an appropriate **MPI** library. The software packages necessary to build the source code of DSP are listed below.

*NOTE: We have also tested installing and running DSP on MacOSX (10.9 or later). This requires the user to install the required packages using Macport. We will prepare the documentation for that soon. For now, please email and consult us for how to install and run DSP on MacOSX.*

### Build essentials

You can install the build essential packages by using the shell script get.essentials in the DSP project repository or by manually installing the packages.

#### Using ./get.essentials script to install the build essential packages

You can run the apt-get commands to install the packages required for DSP by
```bash
sudo ./get.essentials
```

#### Manual installation of the build essential packages

If you have used ./get.essentials, then you can ignore this section and directly go to [External software packages](#External). Before using apt-get, please update the package list by typing
```bash
sudo apt-get update
```
* MPICH -- A version of MPICH is **required** to build and run DSP. It is required to install a version of MPICH, which is available from http://www.mpich.org/downloads/. On Linux, you can also do
```bash
    sudo apt-get install libmpich-dev
```
* CMake -- This package is **required** to build DSP and all the other external software packages. It is available from http://www.cmake.org/download/. On Linux, you can also do
```bash
    sudo apt-get install cmake
```
* BLAS/LAPACK -- These linear algebra libraries are **required** and may already be available on your machine. If not available, the source codes are available from http://www.netlib.org/blas/blas.tgz and http://www.netlib.org/lapack/lapack-3.5.0.tgz. Or on linux you can install them by typing
```bash
    sudo apt-get install libblas-dev
    sudo apt-get install liblapack-dev
```
* svn -- A subversion software package is **required** to download an external package Coin-SMI used in DSP. You can install it by
```bash
    sudo apt-get install subversion
```
* GNU Make/Autoconf/Automake -- These are **required** and likely available on any machine. If not available, you can get them as follows.
```bash
    sudo apt-get install build-essential
    sudo apt-get install autoconf
    sudo apt-get install automake
```
* bzip2/zlib/xtuils-dev -- These packages are **required** to build one of the external software packages used in DSP.
```bash
    sudo apt-get install libbz2-dev
    sudo apt-get install zlib1g-dev
    sudo apt-get install xutils-dev
```

<a name="External"></a>
### External software packages

The following packages are also **required** to build and run DSP and need to be located on ./extra directory before DSP may be built. **DSP will automatically configure and build the external packages once they are located in the right place.**

* MA27 (./extra/ma27-1.0.0) -- MA27 is a library for solving sparse symmetric indefinite linear systems. To build OOQP solver, you must have this installed. This is part of HSL (formerly the Harwell Subroutine Library), a collection of ISO Fortran codes for large scale scientific computation.
  1. Download: http://www.hsl.rl.ac.uk/download/MA27/1.0.0/a/
  2. Unpack the downloaded file.
  3. Move and rename the downloaded directory to ./extras/ma27-1.0.0/

* SCIP Optimization Suite (./extra/scipoptsuite-3.1.1) -- SCIP contains non-commercial solvers for mixed integer programming (MIP) and mixed-integer nonlinear programming (MINLP).
  1. Download: http://scip.zib.de/download.php?fname=scipoptsuite-3.1.1.tgz
  2. Unpack the downloaded file.
  3. Move and rename the downloaded directory to ./extras/scipoptsuite-3.1.1

## Installation

If you have all the prerequisite packages installed on your system, then you need to go to the root directory of DSP and type
```bash
cmake .
```
to configure OOQP. If you wish to install the package in a more permanent location, you may then type
```bash
make install
```
External packages (MA27, OOQP, SCIP Optimization Suite, Smi) used in DSP are built automatically. A shared object is installed in ./lib directory. Once the installation has been successfully done, you need to set environment variable (DY)LD_LIBRARY_PATH.
Please add the following line by replacing \<DSP_SRC_PATH\> with your DSP source directory in ~/.bash_profile (or ~/.bash_aliases): for Linux,
```bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<DSP_SRC_PATH>/lib
```
for Mac,
```bash
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:<DSP_SRC_PATH>/lib
```

## Julia Interface (Required)

DSP uses Julia as a modeling interface for the computational experiments. Julia is a high level dynamic programming language for technical computing, with syntax that is familiar to users of other technical computing environments such as MATLAB and Python. Julia can be downloaded from http://julialang.org/downloads/ If Julia is successfully installed on your machine, then you can start the Julia command-line tool by typing
```bash
julia
```
Now we need to install four Julia packages required to run DSP in the Julia environment. The packages should be installed in the Julia command-line tool. Please update the package list of Julia by typing
```julia
Pkg.update();
```
* [MPI.jl](https://github.com/JuliaParallel/MPI.jl) is an MPI interface package, which can be installed by the following Julia command. MPICH is required for this package.
```julia
    Pkg.add("MPI");
```
* [JuMP.jl](https://github.com/JuliaOpt/JuMP.jl) is a algebraic modeling package in Julia for mathematical programming, which can be installed by the following Julia command:
```julia
    Pkg.add("JuMP");
```
* [StochJuMP.jl](https://github.com/kibaekkim/StochJuMP.jl.git) is a algebraic modeling package in Julia for stochastic programming, which can be installed by the following Julia command:
```julia
    Pkg.clone("https://github.com/kibaekkim/StochJuMP.jl.git");
```
* The [DSPsolver.jl](https://github.com/kibaekkim/DSPsolver.jl.git) package provides an interface to StochJuMP. DSPsolver.jl can be installed by the Julia command
```julia
    Pkg.clone("https://github.com/kibaekkim/DSPsolver.jl.git");
```

##Example (farmer)

An example is provided in subdirectory examples/farmer. You can run the farmer example with the dual decomposition method by typing the command line
```bash
mpirun -n 1 julia examples/farmer/farmer_run.jl DD # dual decomposition
```
You can also solve the problem by using the Benders decomposition and the extensive form solution implemented in DSP:
```bash
julia examples/farmer/farmer_run.jl BD # Benders decomposition
julia examples/farmer/farmer_run.jl DE # Extensive form solution
```

## Reproducing the Computational Results

We provide scripts to reproduce the computational results reported in [the manuscript](http://www.optimization-online.org/DB_FILE/2015/06/4960.pdf). The scripts are available in subdirectory ./experiments/ that contains the bash-shell scripts for running DSP and the MATLAB scripts for reproducing tables and figures. The results will be stored in subdirectory ./output in CSV format.
* Table2.pbs.sh (Table2.sh) -- This runs DSP on cluster with the PBS scheduler for generating Table 2 in the manuscript. If PBS is not available, then Table2.sh can be used on local machine.
 * Table2.m -- This MATLAB script will generate Table2.tex from the CSV files generated by Table2.pbs.sh (or Table2.sh).
 * Figures2and3.m -- This MATLAB script will generate Figure2.eps and Figure3.eps from the CSV files generated by Table2.pbs.sh (or Table2.sh).
 * Figure4.m -- This MATLAB script will generate Figure4.eps from the CSV files generated by Table2.pbs.sh (or Table2.sh).

Similarly, the following scripts run DSP on cluster with the PBS scheduler for generating the figures and the tables in the manuscript.
* Table4.pbs.sh (or Table4.sh) and Table4.m
* Table5.pbs.sh (or Table5.sh) and Table5.m
* Figure5.pbs.sh (or Figure5.sh) and Figure5.m
* Figure6.pbs.sh (or Figure6.sh) and Figure6.m
* Figure8.pbs.sh (or Figure8.sh) and Figure8.m

## Contents of the Distribution

The DSP distribution contains the following top-level subdirectories:
* experiments/ -- All the scripts to reproduce computational results including tables and figures in the manuscript.
* extra/ -- External packages required and used for DSP. These include MA27, OOQP, SCIP and SMI.
* examples/farmer/ -- This contains an example for using DSP with Julia.
* lib/ -- This contains DSP shared object precompiled on a 64-bit linux machine.
* src/ -- Source files for the DSP distribution.

## Credits

DSP has been developed and is maintained by:
* [Kibaek Kim](http://mcs.anl.gov/~kibaekkim/), Mathematics and Computer Science Division, Argonne National Laboratory.
* [Victor M. Zavala](http://mcs.anl.gov/~vzavala/), Mathematics and Computer Science Division, Argonne National Laboratory.

## Publication
* Kibaek Kim and Victor M. Zavala. "[Algorithmic innovations and software for the dual decomposition method applied to stochastic mixed-integer programs](http://www.optimization-online.org/DB_FILE/2015/06/4960.pdf)" Optimization Online, 2015

## Acknowledgements

This material is based upon work supported by the U.S. Department of Energy, Office of Science, under contract number DE-AC02-06CH11357. We gratefully acknowledge the computing resources provided on *Blues*, a high-performance computing cluster operated by the Laboratory Computing Resource Center at Argonne National Laboratory. We thank E. Michael Gertz and Stephen Wright for providing the [OOQP](http://pages.cs.wisc.edu/~swright/ooqp/) software package.
