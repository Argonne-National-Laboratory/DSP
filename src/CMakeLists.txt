cmake_minimum_required(VERSION 2.8)

##################################################
# Set directories
##################################################

set(DEPEND_DIR     "/Users/kibaekkim/Documents/REPOS/DspCoinDepend/build")
set(CPLEX_DIR      "/Users/kibaekkim/Applications/IBM/ILOG/CPLEX_Studio1261/cplex")
set(CPLEX_LIB_DIR  "${CPLEX_DIR}/lib/x86-64_osx/static_pic")
set(CPLEX_INC_DIR  "${CPLEX_DIR}/include/ilcplex")
set(GUROBI_DIR     "/Library/gurobi652")
set(GUROBI_LIB_DIR "${GUROBI_DIR}/mac64/lib")
set(GUROBI_INC_DIR "${GUROBI_DIR}/mac64/include")

##################################################
# add include directories
##################################################

include_directories(
	Model 
	Solver 
	Utility
	${DEPEND_DIR}/include/coin
	${DEPEND_DIR}/include/ooqp
	${CPLEX_INC_DIR}
	${GUROBI_INC_DIR}
)

##################################################
# add files to Dsp library
##################################################

set(DSP_SRC_FILES
	Model/BlkModel.cpp
	Model/DecBlkModel.cpp
	Model/DecTssModel.cpp
	Model/DetModel.cpp
	Model/StoModel.cpp
	Model/TssModel.cpp
	Solver/DantzigWolfe/DwModel.cpp
	Solver/DantzigWolfe/DwMaster.cpp
	Solver/DantzigWolfe/DwBundleDual.cpp
	Solver/DantzigWolfe/DwHeuristic.cpp
	Solver/DantzigWolfe/DwSolverSerial.cpp
	Solver/DantzigWolfe/DwWorker.cpp
	SolverInterface/OoqpEps.cpp
	SolverInterface/OoqpStatus.cpp
	SolverInterface/OsiOoqpSolverInterface.cpp
	TreeSearch/DspModel.cpp
	TreeSearch/DspTreeNode.cpp
	Utility/DspParams.cpp)
if (MPI_CXX_FOUND)
	set(DSP_SRC_FILES ${DSP_SRC_FILES}
		Solver/DantzigWolfe/DwSolverMpi.cpp
		Solver/DantzigWolfe/DwWorkerMpi.cpp
		Utility/DspMpi.cpp)
endif()
add_library(Dsp SHARED
	DspApiEnv.cpp
	DspCInterface.cpp
	${DSP_SRC_FILES})

find_library(SMILIB Smi PATHS ${DEPEND_DIR}/lib)
find_library(OSILIB Osi PATHS ${DEPEND_DIR}/lib)
find_library(COINUTILSLIB CoinUtils PATHS ${DEPEND_DIR}/lib)
find_library(ALPSLIB Alps PATHS ${DEPEND_DIR}/lib)
find_library(OOQPLIB1 ooqpmehrotra PATHS ${DEPEND_DIR}/lib)
find_library(OOQPLIB2 ooqpgendense PATHS ${DEPEND_DIR}/lib)
find_library(OOQPLIB3 ooqpgensparse PATHS ${DEPEND_DIR}/lib)
find_library(OOQPLIB4 ooqpdense PATHS ${DEPEND_DIR}/lib)
find_library(OOQPLIB5 ooqpsparse PATHS ${DEPEND_DIR}/lib)
find_library(OOQPLIB6 ooqpbase PATHS ${DEPEND_DIR}/lib)
find_library(MA27LIB ma27 PATHS /usr/local/lib)
find_library(OSICPXLIB OsiCpx PATHS ${DEPEND_DIR}/lib)
find_library(OSIGRBLIB OsiGrb PATHS ${DEPEND_DIR}/lib)
find_library(CPXLIB cplex PATHS ${CPLEX_LIB_DIR})
find_library(GRBLIB gurobi65 PATHS ${GUROBI_LIB_DIR})

# set DSP libraries
set(DSP_LIBS
	${SMILIB}
	${OSICPXLIB}
	${OSIGRBLIB}
	${OSILIB}
	${COINUTILSLIB}
	${ALPSLIB}
	${CPXLIB}
	${GRBLIB}
	${OOQPLIB1}
	${OOQPLIB2}
	${OOQPLIB3}
	${OOQPLIB4}
	${OOQPLIB5}
	${OOQPLIB6}
	${MA27LIB}
	gfortran
	${LAPACK_LIBRARIES}
)
if (MPI_CXX_FOUND)
	set(DSP_LIBS ${DSP_LIBS} ${MPI_CXX_LIBRARIES})
endif()
target_link_libraries(Dsp ${DSP_LIBS})


##################################################
# INSTALL
##################################################
install(TARGETS Dsp LIBRARY DESTINATION $ENV{PWD}/lib)

##################################################
# SUBDIRECTORY
##################################################
# if (PYDSP_COMPILE)
# 	add_subdirectory(PyDsp)
# endif()

##################################################
# Messages
##################################################
install(CODE "message(\"\")")
install(CODE "message(\"${Blue}** Please add the following line to the .bash_profile or .bashrc in your home directory:${ColourReset}\")")
if(APPLE)
	install(CODE "message(\"${Red}    export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$ENV{PWD}/lib${ColourReset}\")")
else()
	install(CODE "message(\"${Red}    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ENV{PWD}/lib${ColourReset}\")")
endif()
