## Prerequisites

DSP has been tested on several linux distributions and MacOSX (10.9 or later). The software packages necessary to build the source code of DSP are listed below.

### Build Essentials

Installation of DSP requires some libraries. Installing the libraries depends on your distribution and more specifically package managers.

#### Quick Instructions

We provide scripts for the most common distributions. If you cannot find your distribution below, please refer the manual instruction section.

* MacOSX: This requires to have [brew](http://brew.sh) installed.
```bash
sudo ./get.essentials.macosx
```
* Debian based distribution (Debian, Ubuntu, Mint, ...):
```bash
sudo ./get.essentials.debian
```
* Redhat based distribution (Redhat, Fedora, CentOS, ...):
```bash
sudo ./get.essentials.redhat
```
* OpenSUSE:
```bash
sudo ./get.essentials.opensuse
```
* ArchLinux:
```bash
sudo ./get.essentials.archlinux
```

#### Manual Installation

* **CMake** -- This package provides a cross-platform make system. This is available in http://www.cmake.org/download/.
* **BLAS/LAPACK** -- These linear algebra libraries are *required* and may already be available on your machine. If not available, the source codes are available in http://www.netlib.org/blas/blas.tgz and http://www.netlib.org/lapack/lapack-3.5.0.tgz.
* **Subversion** -- A subversion software package is *required* to download an external package Coin-SMI used in DSP. If this is not available, please download the binary from https://subversion.apache.org/packages.html.
* **GNU Make/Autoconf/Automake** -- These are *required* and likely available on any machine. If not available, please get them from the following websites:
  * https://www.gnu.org/software/make/
  * https://www.gnu.org/software/autoconf/autoconf.html
  * https://www.gnu.org/software/automake/
* **bzip2 library** -- This is available in http://www.bzip.org/downloads.html.
* **zlib library** -- This is available in http://www.zlib.net.

## MPI Library ##

MPI library is optional to build and run DSP, but required to run DSP in parallel. To run in parallel, you need to install one of the following libraries.

* **MPICH** -- This is available from http://www.mpich.org.
* **OpenMPI** -- This is available from https://www.open-mpi.org.

## Installation

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

### Binary file

A binary file ``runDsp`` is installed in ./bin directory.

### Shared Object

A shared object is installed in ./lib directory. Once the installation has been successfully done, you need to set environment variable (DY)LD_LIBRARY_PATH.
Please add the following line by replacing \<DSP_SRC_PATH\> with your DSP source directory in ~/.bash_profile (or ~/.bash_aliases): for Linux,
```bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<DSP_SRC_PATH>/lib
```
for Mac,
```bash
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:<DSP_SRC_PATH>/lib
```
