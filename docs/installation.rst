Installation
------------

This guide describes the prerequisites for installation and how to install DSP and Julia interface.

Download
^^^^^^^^

You can clone this repository in your preferred directory by typing::

   git clone https://github.com/kibaekkim/DSP.git

Prerequisites
^^^^^^^^^^^^^

We recommend to install and run DSP on a **Linux** machine with an appropriate **MPI** library. The software packages necessary to build the source code of DSP are listed below.

*NOTE: We have also tested installing and running DSP on MacOSX (10.9 or later). This requires the user to install the required packages using Macport. We will prepare the documentation for that soon. For now, please email and consult us for how to install and run DSP on MacOSX.*

Build Essentials
################

You can install the build essential packages by using the shell script ``./get.essentials`` in the DSP project repository or by manually installing the packages.

Using ./get.essentials script
*****************************

You can run the apt-get commands to install the packages required for DSP::

   sudo ./get.essentials

Manual installation
*******************

If you have used ``./get.essentials``, then you can ignore this section and directly go to `External software packages`_. Before using apt-get, please update the package list::

   sudo apt-get update

* MPICH -- A version of MPICH is **required** to build and run DSP. It is required to install a version of MPICH, which is available from http://www.mpich.org/downloads/. On Linux, you can also do::

    sudo apt-get install libmpich-dev

* CMake -- This package is **required** to build DSP and all the other external software packages. It is available from http://www.cmake.org/download/. On Linux, you can also do::

    sudo apt-get install cmake

* BLAS/LAPACK -- These linear algebra libraries are **required** and may already be available on your machine. If not available, the source codes are available from http://www.netlib.org/blas/blas.tgz and http://www.netlib.org/lapack/lapack-3.5.0.tgz. Or on linux you can install them by typing::

    sudo apt-get install libblas-dev
    sudo apt-get install liblapack-dev

* svn -- A subversion software package is **required** to download an external package Coin-SMI used in DSP. You can install it by::

    sudo apt-get install subversion

* GNU Make/Autoconf/Automake -- These are **required** and likely available on any machine. If not available, you can get them as follows::

    sudo apt-get install build-essential
    sudo apt-get install autoconf
    sudo apt-get install automake

* bzip2/zlib/xtuils-dev -- These packages are **required** to build one of the external software packages used in DSP::

    sudo apt-get install libbz2-dev
    sudo apt-get install zlib1g-dev
    sudo apt-get install xutils-dev

External software packages
##########################

The following packages are also **required** to build and run DSP and need to be located on ./extra directory before DSP may be built. **DSP will automatically configure and build the external packages once they are located in the right place.**

* MA27 -- This is a library for solving sparse symmetric indefinite linear systems. To build OOQP solver, you must have this installed. This is part of HSL (formerly the Harwell Subroutine Library), a collection of ISO Fortran codes for large scale scientific computation.

  1. Download: http://www.hsl.rl.ac.uk/download/MA27/1.0.0/a/
  2. Unpack the downloaded file.
  3. Move and rename the downloaded directory to ``./extras/ma27-1.0.0/``

* SCIP Optimization Suite -- This contains non-commercial solvers for mixed integer programming (MIP) and mixed-integer nonlinear programming (MINLP).

  1. Download: http://scip.zib.de/download.php?fname=scipoptsuite-3.1.1.tgz
  2. Unpack the downloaded file.
  3. Move and rename the downloaded directory to ``./extras/scipoptsuite-3.1.1``

Build and Setting
^^^^^^^^^^^^^^^^^

If you have all the prerequisite packages installed on your system, then you need to go to the root directory of DSP and type::

   cmake .

to configure OOQP. If you wish to install the package in a more permanent location, you may then type::

   make install

External packages (MA27, OOQP, SCIP Optimization Suite, Smi) used in DSP are built automatically. A shared object is installed in ./lib directory. Once the installation has been successfully done, you need to set environment variable ``LD_LIBRARY_PATH``.
Please add the following line by replacing \<DSP_SRC_PATH\> with your DSP source directory in ``~/.bash_profile`` (or ``~/.bash_aliases``)::

   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<DSP_SRC_PATH>/lib

Julia Interface
^^^^^^^^^^^^^^^

DSP uses Julia as a modeling interface for the computational experiments. Julia is a high level dynamic programming language for technical computing, with syntax that is familiar to users of other technical computing environments such as MATLAB and Python. Julia can be downloaded from http://julialang.org/downloads/ If Julia is successfully installed on your machine, then you can start the Julia command-line tool by typing::

   julia

Now we need to install four Julia packages required to run DSP in the Julia environment. The packages should be installed in the Julia command-line tool. Please update the package list of Julia by typing::

   julia> Pkg.update();

The following packages should be installed.

* `MPI.jl <https://github.com/JuliaParallel/MPI.jl>`_ is an MPI interface package, which can be installed by the following Julia command. MPICH is required for this package::

   julia> Pkg.add("MPI");

* `JuMP.jl <https://github.com/JuliaOpt/JuMP.jl>`_ is an algebraic modeling package in Julia for mathematical programming, which can be installed by the following Julia command::

   julia> Pkg.add("JuMP");

* `StochJuMP.jl <https://github.com/kibaekkim/StochJuMP.jl.git>`_ is an algebraic modeling package in Julia for stochastic programming, which can be installed by the following Julia command::

   julia> Pkg.clone("https://github.com/kibaekkim/StochJuMP.jl.git");

* The `DSPsolver.jl <https://github.com/kibaekkim/DSPsolver.jl.git>`_ package provides an interface to StochJuMP. DSPsolver.jl can be installed by the Julia command::

   julia> Pkg.clone("https://github.com/kibaekkim/DSPsolver.jl.git");

