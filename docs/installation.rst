Installation
------------

This describes the following:

* how to correctly clone the repository with submodules
* the prerequisites for installation
* how to install DSP and Julia interface

Download
^^^^^^^^

You can clone this repository in your preferred directory by typing::

   git clone --recursive https://github.com/Argonne-National-Laboratory/DSP.git

or::

   git clone https://github.com/Argonne-National-Laboratory/DSP.git
   cd DSP
   git submodule update --init --recursive


Build and Setting
^^^^^^^^^^^^^^^^^

Please set UserConfig.cmake as follows.

* MA27LIB_DIR is *required* to use OOQP (interior point solver). If you are not using CPLEX, it is recommended to use OOQP for better performance. You can request the library here: http://www.hsl.rl.ac.uk/download/MA27/1.0.0/a/
* It is required to set paths for either CPLEX or SCIP.
  * CPLEX_LIB_DIR is the path to the directory that contains libcplex.a. CPLEX_INC_DIR is the path to the directory that contains cplex.h.
  * If you use SCIP, you need to compile it as a shared library. SCIPOPT_INC_DIR is the path to the SCIP directory. SCIPOPT_LIB_DIR is the path to the shared library (usually ${SCIPOPT_INC_DIR}/lib). SPX_DIR is the path to the SOPLEX directory.
* Once you are set the variables above, please set USER_SETTINGS to ON.

If you have all the prerequisite packages installed on your system, then you need to go to the root directory of DSP and type::

   mkdir build
   cd build
   cmake ..

A shared object is installed in ``./lib`` directory. Once the installation has been successfully done, you need to set environment variable ``(DY)LD_LIBRARY_PATH``.
For Linux::

   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<DSP_SRC_PATH>/lib

Outputs
^^^^^^^

Binary file
***********

A binary file ``runDsp`` is installed in ``./bin`` directory.

Shared object
*************

A shared object is installed in ``./lib`` directory. Once the installation has been successfully done, you need to set environment variable ``(DY)LD_LIBRARY_PATH``.
For Linux::

   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<DSP_SRC_PATH>/lib

For Mac::

   export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:<DSP_SRC_PATH>/lib
