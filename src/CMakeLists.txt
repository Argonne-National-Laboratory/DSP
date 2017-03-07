##################################################
# set directories
##################################################

set(COIN_LIB_DIR "${DEPEND_DIR}/lib")
set(COIN_INC_DIR "${DEPEND_DIR}/include/coin")
set(OOQP_LIB_DIR "${DEPEND_DIR}/lib")
set(OOQP_INC_DIR "${DEPEND_DIR}/include/ooqp")

##################################################
# add include directories
##################################################

include_directories(
	Model 
	Solver 
	Utility
	${COIN_INC_DIR}
	${OOQP_INC_DIR}
	${CPLEX_INC_DIR}
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
add_dependencies(Dsp Alps OOQP)

# set DSP libraries
set(DSP_LIBS
	${COIN_LIB_DIR}/libSmi.a
	${COIN_LIB_DIR}/libOsiCpx.a
	${COIN_LIB_DIR}/libOsi.a
	${COIN_LIB_DIR}/libCoinUtils.a
	${COIN_LIB_DIR}/libAlps.a
	${CPLEX_LIB_DIR}/libcplex.a
	${OOQP_LIB_DIR}/libooqpmehrotra.a
	${OOQP_LIB_DIR}/libooqpgendense.a
	${OOQP_LIB_DIR}/libooqpgensparse.a
	${OOQP_LIB_DIR}/libooqpdense.a
	${OOQP_LIB_DIR}/libooqpsparse.a
	${OOQP_LIB_DIR}/libooqpbase.a
	${MA27LIB_DIR}/libma27.a
	gfortran
	${LAPACK_LIBRARIES}
	${ZLIB_LIBRARIES}
	${BZIP2_LIBRARIES}
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
