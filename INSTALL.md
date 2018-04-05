# Installation

Please set `UserConfig.cmake` as follows.

* `MA27LIB_DIR` is optional to use OOQP (interior point solver). If you are not using CPLEX, it is recommended to use OOQP for better performance. You can request the library here: http://www.hsl.rl.ac.uk/download/MA27/1.0.0/a/
* It is required to set paths for either CPLEX or SCIP. 
  * `CPLEX_LIB_DIR` is the path to the directory that contains `libcplex.a`. `CPLEX_INC_DIR` is the path to the directory that contains `cplex.h`.
  * If you use SCIP, you need to compile it as a shared library. `SCIP_DIR` is the path to the SCIP directory. `SCIP_LIB_DIR`  is the path to the shared library (usually `${SCIP_DIR}/lib`). `SPX_DIR` is the path to the SOPLEX directory.
* Once you are set the variables above, please set `USER_SETTINGS` to `ON`.

If you have all the prerequisite packages installed on your system, then you need to go to the root directory of DSP and type
```bash
mkdir build
cd build
cmake ..
```
to configure OOQP. If you wish to install the package in a more permanent location, you may then type
```bash
make install
```
External packages (MA27, OOQP, SCIP Optimization Suite, Smi) used in DSP are built automatically. 

## Binary file

A binary file ``runDsp`` is installed in ./bin directory.

## Shared Object

A shared object is installed in ./lib directory. Once the installation has been successfully done, you need to set environment variable (DY)LD_LIBRARY_PATH.
Please add the following line by replacing \<DSP_SRC_PATH\> with your DSP source directory in ~/.bash_profile (or ~/.bash_aliases): for Linux,
```bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<DSP_SRC_PATH>/lib
```
for Mac,
```bash
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:<DSP_SRC_PATH>/lib
```
