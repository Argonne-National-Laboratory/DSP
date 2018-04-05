# Prerequisites

DSP has been tested on MacOSX (10.9 or later) and Linux machines. The software packages necessary to build the source code of DSP are listed below.

## Build Essentials

There are some libraries required to run DSP.

### Quick Instructions

If ``apt-get`` is available on your system, please run the following command to install the packages required for DSP. Otherwise, please refer the manual installation.
```bash
sudo ./get.essentials
```

### Manual Installation

* **CMake** -- This package provides a cross-platform make system. This is available in http://www.cmake.org/download/.
* **BLAS/LAPACK** -- These linear algebra libraries are *required* and may already be available on your machine. If not available, the source codes are available in http://www.netlib.org/blas/blas.tgz and http://www.netlib.org/lapack/lapack-3.5.0.tgz.
* **GNU Make/Autoconf/Automake** -- These are *required* and likely available on any machine. If not available, please get them from the following websites:
  * https://www.gnu.org/software/make/
  * https://www.gnu.org/software/autoconf/autoconf.html
  * https://www.gnu.org/software/automake/
* **bzip2 library** -- This is available in http://www.bzip.org/downloads.html.
* **zlib library** -- This is available in http://www.zlib.net.

## External software packages

The following packages are also *required* to build and run DSP and need to be located on ./extra directory before DSP may be built. *DSP will automatically configure and build the external packages once they are located in the right place.*

* **MA27** (./extra/ma27-1.0.0) -- MA27 is a library for solving sparse symmetric indefinite linear systems. To build OOQP solver, you must have this installed. This is part of HSL (formerly the Harwell Subroutine Library), a collection of ISO Fortran codes for large scale scientific computation.
  1. Download: http://www.hsl.rl.ac.uk/download/MA27/1.0.0/a/
  2. Unpack the downloaded file.
  3. Move and rename the downloaded directory to ./extras/ma27-1.0.0/

* **SCIP Optimization Suite** (./extra/scipoptsuite-3.1.1) -- SCIP contains non-commercial solvers for mixed integer programming (MIP) and mixed-integer nonlinear programming (MINLP).
  1. Download: http://scip.zib.de/download.php?fname=scipoptsuite-3.1.1.tgz
  2. Unpack the downloaded file.
  3. Move and rename the downloaded directory to ./extras/scipoptsuite-3.1.1

## MPI Library ##

MPI library is optional to build and run DSP, but required to run DSP in parallel. To run in parallel, you need to install one of the following libraries.

* **MPICH** -- This is available from http://www.mpich.org.
* **OpenMPI** -- This is available from https://www.open-mpi.org.
