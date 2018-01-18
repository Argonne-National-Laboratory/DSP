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
find_package(LAPACK REQUIRED)
find_package(ZLIB REQUIRED)
find_package(BZip2 REQUIRED)
find_package(MPI)
find_package(OpenMP)

##################################################
# My include files
##################################################
include(UserConfig.cmake)
include(color_scheme.cmake)

##################################################
# Compiler settings
##################################################

# if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
# 	message(FATAL_ERROR "We currently do not support Clang compiler")
# elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
# 	message(FATAL_ERROR "We currently do not support Intel compiler")
# elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
# 	message(FATAL_ERROR "We currently do not support MSVC compiler")
# endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++0x -DNDEBUG")
if(MPI_CXX_FOUND)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS} -DDSP_HAS_MPI")
endif()
if (OPENMP_FOUND)
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
if(APPLE)
	set(CMAKE_MACOSX_RPATH OFF)
endif()

##################################################
# add include directories
##################################################
include_directories(src)
if (MPI_CXX_FOUND)
	include_directories(${MPI_CXX_INCLUDE_PATH})
endif()

##################################################
# SUBDIRECTORY
##################################################
set(DEPEND_DIR $ENV{PWD})
add_subdirectory(extra)
add_subdirectory(src)
