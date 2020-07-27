cmake_minimum_required(VERSION 3.3)

project(farmers C)

option(CXX "enable C++ compilation" ON)
if(CXX)
    enable_language(CXX)
endif()

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Debug")
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# list source files here
set(sources farmers.cpp)

add_executable(farmers ${sources})

# include files DspCInterface.h and all files in /Users/yingqiuzhang/Documents/DSPtest/DSP/build/include/coin
include_directories(${CMAKE_PROJECT_NAME} " ")
include_directories(${CMAKE_PROJECT_NAME} " ")
# link the libDsp.dylib
target_link_libraries(${CMAKE_PROJECT_NAME} " ")