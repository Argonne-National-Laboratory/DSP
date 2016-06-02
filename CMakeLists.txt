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

# Check gcc version
if (CMAKE_COMPILER_IS_GNUCC)
	execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
	string(REGEX MATCHALL "[0-9]+" GCC_VERSION_COMPONENTS ${GCC_VERSION})
	list(GET GCC_VERSION_COMPONENTS 0 GCC_MAJOR)
	list(GET GCC_VERSION_COMPONENTS 1 GCC_MINOR)
endif()

if (DSP_ASSERT EQUAL 1)
	set(DSP_NDEBUG "")
else()
	set(DSP_NDEBUG "-DNDEBUG")
endif()

if (${GCC_MAJOR} LESS 5 AND ${GCC_MINOR} LESS 7)
	set(CMAKE_CXX_FLAGS "${MPI_CXX_COMPILE_FLAGS} -O3 -fopenmp -std=c++0x ${DSP_NDEBUG} -DUSE_OMP")
else()
	set(CMAKE_CXX_FLAGS "${MPI_CXX_COMPILE_FLAGS} -O3 -fopenmp -std=c++11 ${DSP_NDEBUG} -DUSE_OMP")
endif()

if(APPLE)
	set(CMAKE_MACOSX_RPATH OFF)
endif()

# Set path variables
set(EXTRA_DIR ${CMAKE_SOURCE_DIR}/extra)

# Build and install external packages
include(ExternalProject)

if (MA27_DIR)
	message(STATUS "MA27_DIR is set at ${MA27_DIR}")
else()
	set(MA27_DIR ${EXTRA_DIR}/ma27-1.0.0)
	# ma27-1.0.0
	ExternalProject_Add(ma27-1.0.0
		PREFIX ${MA27_DIR}
		SOURCE_DIR ${MA27_DIR}
		CONFIGURE_COMMAND ${MA27_DIR}/configure --prefix=${MA27_DIR} FFLAGS=-fPIC
		BINARY_DIR ${MA27_DIR}
		BUILD_COMMAND make
		INSTALL_COMMAND make install
	)
endif()

if (OOQP_DIR)
	message(STATUS "OOQP_DIR is set at ${OOQP_DIR}")
else()
	set(OOQP_DIR ${EXTRA_DIR}/OOQP-0.99.25)
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
endif()

if (SCIPOPTSUITE_DIR)
	message(STATUS "SCIPOPTSUITE_DIR is set at ${SCIPOPTSUITE_DIR}")
	set(SCIP_DIR ${SCIPOPTSUITE_DIR}/scip-3.1.1)
	set(SOPLEX_DIR ${SCIPOPTSUITE_DIR}/soplex-2.0.1)
else()
	set(SCIPOPTSUITE_DIR ${EXTRA_DIR}/scipoptsuite-3.1.1)
	set(SCIP_DIR ${SCIPOPTSUITE_DIR}/scip-3.1.1)
	set(SOPLEX_DIR ${SCIPOPTSUITE_DIR}/soplex-2.0.1)
	# SCIP Optimization Suite (3.1.1)
	ExternalProject_Add(scipoptsuite-3.1.1
		PREFIX ${SCIPOPTSUITE_DIR}
		SOURCE_DIR ${SCIPOPTSUITE_DIR}
		CONFIGURE_COMMAND ""
		BINARY_DIR ${SCIPOPTSUITE_DIR}
		BUILD_COMMAND make scipoptlib SHARED=true GMP=false ZLIB=false READLINE=false
		INSTALL_COMMAND ""
	)
endif()

if (SMI_DIR)
	message(STATUS "SMI_DIR is set at ${SMI_DIR}")
else()
	set(SMI_DIR ${EXTRA_DIR}/Smi)
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
endif()

#set(LAZYCPLEX_INCDIR "/opt/ibm/ILOG/CPLEX_Studio1261/cplex/include/ilcplex")
#set(LAZYCPLEX_LIB "-L/opt/ibm/ILOG/CPLEX_Studio1261/cplex/lib/x86-64_linux/static_pic -lcplex -lm -lpthread")
#set(SMI_FLAGS "-L${SOPLEX_DIR}/lib -lsoplex")
#ExternalProject_Add(Smi
#	DEPENDS scipoptsuite-3.1.1
#	PREFIX ${SMI_DIR}
#	DOWNLOAD_DIR ${EXTRA_DIR}
#	SVN_REPOSITORY https://projects.coin-or.org/svn/Smi/stable/0.95
#	SVN_TRUST_CERT 1
#	UPDATE_COMMAND ""
#	SOURCE_DIR ${SMI_DIR}
#	CONFIGURE_COMMAND ${SMI_DIR}/configure --enable-static --disable-shared --with-cplex-incdir=${LAZYCPLEX_INCDIR} --with-cplex-lib=${LAZYCPLEX_LIB} --with-soplex-incdir=${SOPLEX_DIR}/src --with-soplex-lib=${SMI_FLAGS} ADD_CXXFLAGS=-fPIC CXX=${CMAKE_CXX_COMPILER}
#	BINARY_DIR ${SMI_DIR}
#	BUILD_COMMAND make
#	INSTALL_COMMAND make install
#)

# SMI installation script does not copy SmiQuadratic.hpp. (their bug)
add_custom_command(TARGET Smi POST_BUILD COMMAND cp ${SMI_DIR}/Smi/src/SmiQuadratic.hpp ${SMI_DIR}/include/coin/)

# add include directories
include_directories(
	src
	src/Model
	src/Solver
	src/Utility
	src/SolverInterface
	${LAZYCPLEX_INCDIR}
	${OOQP_DIR}/include
	${SCIP_DIR}/src
	${SCIP_DIR}/src/scip
	${SOPLEX_DIR}/src
	${SMI_DIR}/include/coin
	${MPI_CXX_INCLUDE_PATH}
)

# compile DSP source code
add_subdirectory(src)

