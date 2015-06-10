cmake_minimum_required(VERSION 2.8)

# Should use GNU compilers for unix-type OS, including Mac and CygWin
if(UNIX)
	set(CMAKE_C_COMPILER "gcc")
	set(CMAKE_CXX_COMPILER "g++")
	set(CMAKE_Fortran_COMPILER "gfortran")
endif()

# Project name
project(DSP)
message("The system information: ${CMAKE_SYSTEM}.${CMAKE_SYSTEM_PROCESSOR}")

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
	message(FATAL_ERROR "We currently do not support Clang compiler")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
	message(FATAL_ERROR "We currently do not support Intel compiler")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
	message(FATAL_ERROR "We currently do not support MSVC compiler")
endif()

# REQUIRED PACKAGES
find_package(MPI REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(ZLIB REQUIRED)
find_package(BZip2 REQUIRED)

set(CMAKE_CXX_FLAGS "${MPI_CXX_COMPILE_FLAGS} -O3 -fopenmp -DNDEBUG -DUSE_OMP")

if(APPLE)
	set(CMAKE_MACOSX_RPATH OFF)
endif()

# Set path variables
set(EXTRA_DIR ${CMAKE_SOURCE_DIR}/extra)
set(MA27_DIR ${EXTRA_DIR}/ma27-1.0.0)
set(OOQP_DIR ${EXTRA_DIR}/OOQP-0.99.25)
set(SCIPOPTSUITE_DIR ${EXTRA_DIR}/scipoptsuite-3.1.1)
set(SCIP_DIR ${SCIPOPTSUITE_DIR}/scip-3.1.1)
set(SOPLEX_DIR ${SCIPOPTSUITE_DIR}/soplex-2.0.1)
set(SMI_DIR ${EXTRA_DIR}/Smi)

# Build and install external packages
include(ExternalProject)

# ma27-1.0.0
ExternalProject_Add(ma27-1.0.0
	PREFIX ${MA27_DIR}
	SOURCE_DIR ${MA27_DIR}
	CONFIGURE_COMMAND ${MA27_DIR}/configure --prefix=${MA27_DIR} FFLAGS=-fPIC
	BINARY_DIR ${MA27_DIR}
	BUILD_COMMAND make
	INSTALL_COMMAND make install
)

# OOQP-0.99.25
set(OOQP_BLAS "BLAS=${BLAS_LIBRARIES}")
string(REPLACE ";" "\ " OOQP_BLAS "${OOQP_BLAS}")
ExternalProject_Add(OOQP-0.99.25
	DEPENDS ma27-1.0.0
	PREFIX ${OOQP_DIR}
	SOURCE_DIR ${OOQP_DIR}
	CONFIGURE_COMMAND ${OOQP_DIR}/configure --prefix=${OOQP_DIR} CXXFLAGS=-fPIC FFLAGS=-fPIC MA27LIB=${MA27_DIR}/lib/libma27.a ${OOQP_BLAS} CXX=${CMAKE_CXX_COMPILER}
	BINARY_DIR ${OOQP_DIR}
	BUILD_COMMAND make
	INSTALL_COMMAND make install_headers
)

# SCIP Optimization Suite (3.1.1)
ExternalProject_Add(scipoptsuite-3.1.1
	PREFIX ${SCIPOPTSUITE_DIR}
	SOURCE_DIR ${SCIPOPTSUITE_DIR}
	CONFIGURE_COMMAND ""
	BINARY_DIR ${SCIPOPTSUITE_DIR}
	BUILD_COMMAND make scipoptlib SHARED=true GMP=false ZLIB=false READLINE=false
	INSTALL_COMMAND ""
)

# COIN-OR Smi (0.95)
set(SMI_FLAGS "-L${SOPLEX_DIR}/lib -lsoplex")
ExternalProject_Add(Smi
	DEPENDS scipoptsuite-3.1.1
	PREFIX ${SMI_DIR}
	DOWNLOAD_DIR ${EXTRA_DIR}
	SVN_REPOSITORY https://projects.coin-or.org/svn/Smi/stable/0.95
	SVN_TRUST_CERT 1
	UPDATE_COMMAND ""
	SOURCE_DIR ${SMI_DIR}
	CONFIGURE_COMMAND ${SMI_DIR}/configure --enable-static --disable-shared --with-soplex-incdir=${SOPLEX_DIR}/src --with-soplex-lib=${SMI_FLAGS} ADD_CXXFLAGS=-fPIC CXX=${CMAKE_CXX_COMPILER}
	BINARY_DIR ${SMI_DIR}
	BUILD_COMMAND make
	INSTALL_COMMAND make install
)

# SMI installation script does not copy SmiQuadratic.hpp. (their bug)
add_custom_command(TARGET Smi POST_BUILD COMMAND cp ${SMI_DIR}/Smi/src/SmiQuadratic.hpp ${SMI_DIR}/include/coin/)

# add include directories
include_directories(
	src
	src/Model
	src/Solver
	src/Utility
	${OOQP_DIR}/include
	${SCIP_DIR}/src
	${SCIP_DIR}/src/scip
	${SOPLEX_DIR}/src
	${SMI_DIR}/include/coin
	${MPI_CXX_INCLUDE_PATH}
)

# compile DSP source code
add_subdirectory(src)

