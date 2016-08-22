cmake_minimum_required(VERSION 2.8)

# Should use GNU compilers for unix-type OS, including Mac and CygWin
if(UNIX)
	set(CMAKE_C_COMPILER "gcc")
	set(CMAKE_CXX_COMPILER "g++")
	set(CMAKE_Fortran_COMPILER "gfortran")
endif()

##################################################
# Project name
##################################################
project(DSP)
message(STATUS "The system information: ${CMAKE_SYSTEM}.${CMAKE_SYSTEM_PROCESSOR}")

##################################################
# REQUIRED PACKAGES
##################################################
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(ZLIB REQUIRED)
find_package(BZip2 REQUIRED)
find_package(MPI)
find_package(PythonLibs)

##################################################
# Build and install external packages
##################################################
include(ExternalProject)

##################################################
# My include files
##################################################
include(color_scheme.cmake)
include(config_path.cmake)

##################################################
# Compiler settings
##################################################

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
	message(FATAL_ERROR "We currently do not support Clang compiler")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
	message(FATAL_ERROR "We currently do not support Intel compiler")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
	message(FATAL_ERROR "We currently do not support MSVC compiler")
endif()

# Check gcc version
if (CMAKE_COMPILER_IS_GNUCC)
	execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
	string(REGEX MATCHALL "[0-9]+" GCC_VERSION_COMPONENTS ${GCC_VERSION})
	list(GET GCC_VERSION_COMPONENTS 0 GCC_MAJOR)
	list(GET GCC_VERSION_COMPONENTS 1 GCC_MINOR)
endif()

# Set flags
if (${GCC_MAJOR} LESS 5 AND ${GCC_MINOR} LESS 7)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif()
if(MPI_CXX_FOUND)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS} -O3 -DNDEBUG -DDSP_HAS_MPI")
else()
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG")
endif()
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -DUSE_OMP")

if(APPLE)
	set(CMAKE_MACOSX_RPATH OFF)
endif()

##################################################
# PYTHON INTERFACE
##################################################
set(PYDSP_COMPILE OFF)
if (PYTHONLIBS_FOUND)
	find_package(SWIG)
	if (SWIG_FOUND)
		if (MPI_CXX_FOUND)
			execute_process(COMMAND python ${CMAKE_SOURCE_DIR}/get_mpi4py_include.py
				RESULT_VARIABLE MPI4PY_NOT_FOUND
				OUTPUT_VARIABLE MPI4PY_INCLUDE)
			if (MPI4PY_NOT_FOUND EQUAL 0)
				set(PYDSP_COMPILE ON)
				set(PYDSP_I "${CMAKE_SOURCE_DIR}/src/PyDsp/PyDspMpi.i")
				string(REPLACE "\n" "" MPI4PY_INCLUDE ${MPI4PY_INCLUDE})
			else()
				message(WARNING "${Red}Python interface requires mpi4py python module.${ColourReset}")
			endif()
		else()
			set(PYDSP_COMPILE ON)
			set(PYDSP_I "${CMAKE_SOURCE_DIR}/src/PyDsp/PyDsp.i")
		endif()
	else()
		message(WARNING "${Red}Python interface requires SWIG.${ColourReset}")
	endif()
else()
	message(WARNING "${Red}Python interface requires phyton....${ColourReset}")
endif()
if (PYDSP_COMPILE)
	message(STATUS "${Blue}Python interface will be compiled.${ColourReset}")
	execute_process(COMMAND mkdir -p $ENV{PWD}/src)
	execute_process(COMMAND mkdir -p $ENV{PWD}/src/PyDsp)
	execute_process(COMMAND swig -python -c++ -I${MPI4PY_INCLUDE} -includeall -o $ENV{PWD}/src/PyDsp/PyDsp.cpp ${PYDSP_I})
else()
	message(WARNING "${Red}Python interface will not be compiled.${ColourReset}")
endif()

##################################################
# option for enabling ma57 (not working well)
##################################################
OPTION(HAS_MA57 "HAS MA57" OFF)

##################################################
# BUILD EXTERNAL PACKAGES
##################################################
set(EXTRA_DIR ${CMAKE_SOURCE_DIR}/extra)
if (READ_CONF)
	message(STATUS "MA27_DIR is set at ${MA27_DIR}")
	message(STATUS "OOQP_DIR is set at ${OOQP_DIR}")
	message(STATUS "SCIPOPTSUITE_DIR is set at ${SCIPOPTSUITE_DIR}")
	set(SCIP_DIR ${SCIPOPTSUITE_DIR}/scip-3.1.1)
	set(SOPLEX_DIR ${SCIPOPTSUITE_DIR}/soplex-2.0.1)
	message(STATUS "SMI_DIR is set at ${SMI_DIR}")
else(READ_CONF)
	# OOQP-0.99.25
	set(OOQP_DIR ${EXTRA_DIR}/OOQP-0.99.25)
	set(OOQP_BLAS "BLAS=${BLAS_LIBRARIES}")
	string(REPLACE ";" "\ " OOQP_BLAS "${OOQP_BLAS}")

	if (HAS_MA57)
		# compile flags
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDSP_HAS_MA57")
		# ma57-3.9.0
		set(MA57_DIR ${EXTRA_DIR}/ma57-3.9.0)
		ExternalProject_Add(ma57
			#DEPENDS metis
			PREFIX ${MA57_DIR}
			SOURCE_DIR ${MA57_DIR}
			CONFIGURE_COMMAND ${MA57_DIR}/configure --prefix=${MA57_DIR} --with-metis=${METIS_DIR}/libmetis.a FFLAGS=-fPIC
			BINARY_DIR ${MA57_DIR}
			BUILD_COMMAND make
			INSTALL_COMMAND make install)
		# metis-4.0.3
		set(METIS_DIR ${EXTRA_DIR}/metis-4.0.3)
		add_custom_command(TARGET ma57 PRE_BUILD
			COMMAND make COPTIONS=-fPIC
			WORKING_DIRECTORY ${METIS_DIR})
		# OOQP
		set(OOQP_DEPENDS "ma57")
		set(OOQP_HSL_LIB "MA57LIB=${MA57_DIR}/lib/libma57.a")
		execute_process(COMMAND cp ${OOQP_DIR}/src/QpGen/Makefile.inc.ma57 ${OOQP_DIR}/src/QpGen/Makefile.inc)
	else(HAS_MA57)
		# ma27-1.0.0
		set(MA27_DIR ${EXTRA_DIR}/ma27-1.0.0)
		ExternalProject_Add(ma27
			PREFIX ${MA27_DIR}
			SOURCE_DIR ${MA27_DIR}
			CONFIGURE_COMMAND ${MA27_DIR}/configure --prefix=${MA27_DIR} FFLAGS=-fPIC
			BINARY_DIR ${MA27_DIR}
			BUILD_COMMAND make
			INSTALL_COMMAND make install)
		# OOQP
		set(OOQP_DEPENDS "ma27")
		set(OOQP_HSL_LIB "MA27LIB=${MA27_DIR}/lib/libma27.a")
		execute_process(COMMAND cp ${OOQP_DIR}/src/QpGen/Makefile.inc.ma27 ${OOQP_DIR}/src/QpGen/Makefile.inc)
	endif()
	# OOQP
	ExternalProject_Add(OOQP-0.99.25
		DEPENDS ${OOQP_DEPENDS}
		PREFIX ${OOQP_DIR}
		SOURCE_DIR ${OOQP_DIR}
		CONFIGURE_COMMAND ${OOQP_DIR}/configure --prefix=${OOQP_DIR} CXXFLAGS=-fPIC FFLAGS=-fPIC ${OOQP_HSL_LIB} ${OOQP_BLAS} CXX=${CMAKE_CXX_COMPILER}
		BINARY_DIR ${OOQP_DIR}
		BUILD_COMMAND make all_libs
		INSTALL_COMMAND make install_headers)

	# SCIP Optimization Suite (3.1.1)
	set(SCIPOPTSUITE_DIR ${EXTRA_DIR}/scipoptsuite-3.1.1)
	set(SCIP_DIR ${SCIPOPTSUITE_DIR}/scip-3.1.1)
	set(SOPLEX_DIR ${SCIPOPTSUITE_DIR}/soplex-2.0.1)
	ExternalProject_Add(scipoptsuite
		PREFIX ${SCIPOPTSUITE_DIR}
		SOURCE_DIR ${SCIPOPTSUITE_DIR}
		CONFIGURE_COMMAND ""
		BINARY_DIR ${SCIPOPTSUITE_DIR}
		BUILD_COMMAND make scipoptlib SHARED=true GMP=false ZLIB=false READLINE=false
		INSTALL_COMMAND "")
	
	# COIN-OR Smi (0.95)
	set(SMI_DIR ${EXTRA_DIR}/Smi)
	set(SMI_FLAGS "-L${SOPLEX_DIR}/lib -lsoplex")
	ExternalProject_Add(Smi
		DEPENDS scipoptsuite
		PREFIX ${SMI_DIR}
		DOWNLOAD_DIR ${EXTRA_DIR}
		SVN_REPOSITORY https://projects.coin-or.org/svn/Smi/stable/0.95
		SVN_TRUST_CERT 1
		UPDATE_COMMAND ""
		SOURCE_DIR ${SMI_DIR}
		CONFIGURE_COMMAND ${SMI_DIR}/configure --enable-static --disable-shared --with-soplex-incdir=${SOPLEX_DIR}/src --with-soplex-lib=${SMI_FLAGS} ADD_CXXFLAGS=-fPIC CXX=${CMAKE_CXX_COMPILER}
		BINARY_DIR ${SMI_DIR}
		BUILD_COMMAND make
		INSTALL_COMMAND make install)
	# SMI installation script does not copy SmiQuadratic.hpp. (their bug)
	add_custom_command(TARGET Smi POST_BUILD COMMAND cp ${SMI_DIR}/Smi/src/SmiQuadratic.hpp ${SMI_DIR}/include/coin/)
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

##################################################
# add include directories
##################################################
set(DSP_INCLUDE_DIRS
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
	${MPI_CXX_INCLUDE_PATH})
if (MPI_CXX_FOUND)
	set(DSP_INCLUDE_DIRS ${DSP_INCLUDE_DIRS} ${MPI_CXX_INCLUDE_PATH})
endif()
if (PYTHONLIBS_FOUND)
	set(DSP_INCLUDE_DIRS ${DSP_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS} ${MPI4PY_INCLUDE})
	if (MPI_CXX_FOUND)
		set(DSP_INCLUDE_DIRS ${DSP_INCLUDE_DIRS} ${MPI4PY_INCLUDE})
	endif()
endif()
include_directories(${DSP_INCLUDE_DIRS})

##################################################
# SUBDIRECTORY
##################################################
add_subdirectory(src)
