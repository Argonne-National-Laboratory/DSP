Installation
------------

This guide describes the prerequisites for installation and how to install DSP and Julia interface.

Download
^^^^^^^^

You can clone this repository in your preferred directory by typing::

   git clone --recursive https://github.com/Argonne-National-Laboratory/DSP.git

Prerequisites
^^^^^^^^^^^^^

We recommend to install and run DSP on a **Linux** machine with an appropriate **MPI** library. The software packages necessary to build the source code of DSP are listed below.

*NOTE: We have also tested installing and running DSP on MacOSX (10.9 or later). We will prepare the documentation for that soon. For now, please email and consult us for how to install and run DSP on MacOSX.*

Build Essentials
################

You can install the build essential packages by using the shell script ``./get.essentials`` in the DSP project repository or by manually installing the packages.

Using ./get.essentials script
*****************************

You can run the apt-get commands to install the packages required for DSP::

   sudo ./get.essentials

Manual installation
*******************

If you have used ``./get.essentials``, then you can ignore this section and directly go to `External software packages`_.

* CMake -- This package is **required** to build DSP and all the other external software packages. It is available from http://www.cmake.org/download/.
* BLAS/LAPACK -- These linear algebra libraries are **required** and may already be available on your machine. If not available, the source codes are available from http://www.netlib.org/blas/blas.tgz and http://www.netlib.org/lapack/lapack-3.5.0.tgz.
* GNU Make/Autoconf/Automake -- These are **required** and likely available on any machine.
* bzip2/zlib/xtuils-dev -- These packages are **required** to build one of the external software packages used in DSP.

MPI Library
***********

MPI library is optional to build and run DSP, but required to run DSP in parallel. To run in parallel, you need to install one of the following libraries.

Build and Setting
^^^^^^^^^^^^^^^^^

Please set UserConfig.cmake as follows.

* MA27LIB_DIR is optional to use OOQP (interior point solver). If you are not using CPLEX, it is recommended to use OOQP for better performance. You can request the library here: http://www.hsl.rl.ac.uk/download/MA27/1.0.0/a/
* It is required to set paths for either CPLEX or SCIP.
  * CPLEX_LIB_DIR is the path to the directory that contains libcplex.a. CPLEX_INC_DIR is the path to the directory that contains cplex.h.
  * If you use SCIP, you need to compile it as a shared library. SCIP_DIR is the path to the SCIP directory. SCIP_LIB_DIR is the path to the shared library (usually ${SCIP_DIR}/lib). SPX_DIR is the path to the SOPLEX directory.
* Once you are set the variables above, please set USER_SETTINGS to ON.

If you have all the prerequisite packages installed on your system, then you need to go to the root directory of DSP and type::

   mkdir build
   cd build
   cmake ..

to configure OOQP. If you wish to install the package in a more permanent location, you may then type::

   make install

Binary file
^^^^^^^^^^^

A binary file ``runDsp`` is installed in ``./bin`` directory.

Shared object
^^^^^^^^^^^^^

A shared object is installed in ``./lib`` directory. Once the installation has been successfully done, you need to set environment variable ``(DY)LD_LIBRARY_PATH``.
For Linux::

   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<DSP_SRC_PATH>/lib

For Mac::

   export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:<DSP_SRC_PATH>/lib

Julia Interface
^^^^^^^^^^^^^^^

DSP uses Julia as a modeling interface for the computational experiments. Julia is a high level dynamic programming language for technical computing, with syntax that is familiar to users of other technical computing environments such as MATLAB and Python. Julia can be downloaded from http://julialang.org/downloads/ If Julia is successfully installed on your machine, then you can start the Julia command-line tool by typing::

   julia

Now we need to install four Julia packages required to run DSP in the Julia environment. The packages should be installed in the Julia command-line tool. Please update the package list of Julia by typing::

   julia> Pkg.update();

The following packages should be installed.

* `JuMP.jl <https://github.com/JuliaOpt/JuMP.jl>`_ is an algebraic modeling package in Julia for mathematical programming, which can be installed by the following Julia command::

   julia> Pkg.add("JuMP");

* The `Dsp.jl <https://github.com/kibaekkim/Dsp.jl.git>`_ package provides an interface to JuMP. Dsp.jl can be installed by the Julia command::

   julia> Pkg.clone("https://github.com/kibaekkim/Dsp.jl.git");

This package is optional.
* `MPI.jl <https://github.com/JuliaParallel/MPI.jl>`_ is an MPI interface package, which can be installed by the following Julia command. MPICH is required for this package::

   julia> Pkg.add("MPI");

