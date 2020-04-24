cmake_minimum_required(VERSION 3.5)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(color_scheme.cmake)
include(UserConfig.cmake)
set(DEPEND_DIR $ENV{PWD})

project(DSP)
message(STATUS "The system information: ${CMAKE_SYSTEM}.${CMAKE_SYSTEM_PROCESSOR}")

enable_language(CXX)
enable_language(Fortran)

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(ZLIB REQUIRED)
find_package(BZip2 REQUIRED)
find_package(MPI)

if (MPI_CXX_FOUND)
	set(CMAKE_C_COMPILER "mpicc")
	set(CMAKE_CXX_COMPILER "mpicxx")
	set(CMAKE_Fortran_COMPILER "mpif90")
endif()

if(APPLE)
	set(LIBRARY_PATH_VAR_NAME "DYLD_LIBRARY_PATH")
else(APPLE)
	set(LIBRARY_PATH_VAR_NAME "LD_LIBRARY_PATH")
endif(APPLE)

find_library(GFORTRANLIB NAMES gfortran libgfortran PATHS ENV ${LIBRARY_PATH_VAR_NAME})
if(NOT GFORTRANLIB)
	message(FATAL_ERROR "gfortran library is required, but cannot be found. Please either locate the library to a searchable path or set a proper environment variable for the path.")
endif(NOT GFORTRANLIB)

find_library(SCIPLIB scip PATHS ${SCIPOPT_LIB_DIR})
find_library(MA27LIB  ma27  PATHS ${MA27LIB_DIR})
find_library(CPLEXLIB cplex PATHS ${CPLEX_LIB_DIR})

# if(APPLE)
# 	set(CMAKE_MACOSX_RPATH OFF)
# endif()

add_subdirectory(extra)
# add_subdirectory(src)

# enable_testing()
# add_subdirectory(test)
