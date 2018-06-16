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
	Model Solver Utility
	${COIN_INC_DIR}
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
	Solver/DantzigWolfe/DwBranchInt.cpp
	Solver/DantzigWolfe/DwBranchNonant.cpp
	Solver/DantzigWolfe/DwBranchNonant2.cpp
	Solver/DantzigWolfe/DwBranchGenDisj.cpp
	Solver/DantzigWolfe/DwModel.cpp
	Solver/DantzigWolfe/DwMaster.cpp
	Solver/DantzigWolfe/DwBundleDual.cpp
	Solver/DantzigWolfe/DwHeuristic.cpp
	Solver/DantzigWolfe/DwSolverSerial.cpp
	Solver/DantzigWolfe/DwWorker.cpp
	#SolverInterface/OoqpEps.cpp
	#SolverInterface/OoqpStatus.cpp
	#SolverInterface/OsiOoqpSolverInterface.cpp
	TreeSearch/DspModel.cpp
	TreeSearch/DspTreeNode.cpp
	Utility/DspParams.cpp)
if (MPI_CXX_FOUND)
	set(DSP_SRC_FILES ${DSP_SRC_FILES}
		Solver/DantzigWolfe/DwSolverMpi.cpp
		Solver/DantzigWolfe/DwWorkerMpi.cpp
		Utility/DspMpi.cpp)
endif()

##################################################
# add link directories
##################################################

link_directories(
	${COIN_LIB_DIR} ${CPLEX_LIB_DIR}
)

##################################################
# add libraries
##################################################
set(DSP_LIBS
	libSmi.a libOsi.a libCoinUtils.a libAlps.a
	${LAPACK_LIBRARIES}
	libOsiCpx.a cplex
	gfortran
	${ZLIB_LIBRARIES}
	${BZIP2_LIBRARIES}
	-ldl
)
if (MPI_CXX_FOUND)
	set(DSP_LIBS ${DSP_LIBS} ${MPI_CXX_LIBRARIES})
endif()

##################################################
# Do we want to use PIPS-IPM?
##################################################

#find_library(PIPS_LIB pipsipm-shared PATHS ${PIPS_LIB_DIR})
#if(${PIPS_LIB} MATCHES "PIPS_LIB-NOTFOUND")
#else()
#	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAS_PIPS")
#	include_directories(
#		${PIPS_DIR}/Input
#		${PIPS_DIR}/PIPS-IPM/Core/Abstract
#		${PIPS_DIR}/PIPS-IPM/Core/Vector
#		${PIPS_DIR}/PIPS-IPM/Core/QpGen
#		${PIPS_DIR}/PIPS-IPM/Core/QpStoch
#		${PIPS_DIR}/PIPS-IPM/Core/QpSolvers
#		${PIPS_DIR}/PIPS-IPM/Core/StochLinearAlgebra
#		${PIPS_DIR}/PIPS-IPM/Core/SparseLinearAlgebra
#	)
#	set(DSP_SRC_FILES ${DSP_SRC_FILES}
#		${PIPS_DIR}/Input/stochasticInput.cpp
#		Solver/DantzigWolfe/DwBundleDualPips.cpp
#		Solver/DantzigWolfe/DwWorkerPips.cpp
#		Solver/DantzigWolfe/PipsInterface.cpp
#	)
#	link_directories(${PIPS_LIB_DIR})
#	set(DSP_LIBS ${DSP_LIBS}
#		pipsipm-shared
#		ma57 metis pardiso
#	)
#endif()

##################################################
# Shared library
##################################################

add_library(DspCore ${DSP_SRC_FILES})
set_property(TARGET DspCore PROPERTY POSITION_INDEPENDENT_CODE ON)
add_library(Dsp SHARED DspApiEnv.cpp DspCInterface.cpp)
add_executable(runDsp dsp.cpp DspApiEnv.cpp DspCInterface.cpp)
add_dependencies(Dsp Alps OOQP)
add_dependencies(runDsp Alps OOQP)

if (APPLE)
	target_link_libraries(Dsp DspCore ${DSP_LIBS})
	target_link_libraries(runDsp DspCore ${DSP_LIBS})
else()
	target_link_libraries(Dsp DspCore ${DSP_LIBS} "-Wl,--no-undefined")
	target_link_libraries(runDsp DspCore ${DSP_LIBS} "-Wl,--no-undefined")
endif(APPLE)

##################################################
# INSTALL
##################################################
install(TARGETS Dsp LIBRARY DESTINATION $ENV{PWD}/lib)
install(TARGETS runDsp DESTINATION $ENV{PWD}/bin)

##################################################
# Messages
##################################################
install(CODE "message(\"\")")
install(CODE "message(\"${Blue}** Please add the following line to the .bash_profile or .bashrc in your home directory:${ColourReset}\")")
set(LIBRARY_PATH_VAR_NAME "LD_LIBRARY_PATH")
if(APPLE)
	set(LIBRARY_PATH_VAR_NAME "DYLD_LIBRARY_PATH")
else()
endif()
if(${PIPS_LIB} MATCHES "PIPS_LIB-NOTFOUND")
	install(CODE "message(\"${Red}    export ${LIBRARY_PATH_VAR_NAME}=$${LIBRARY_PATH_VAR_NAME}:$ENV{PWD}/lib${ColourReset}\")")
else()
	install(CODE "message(\"${Red}    export ${LIBRARY_PATH_VAR_NAME}=$${LIBRARY_PATH_VAR_NAME}:$ENV{PWD}/lib:${PIPS_LIB_DIR}${ColourReset}\")")
endif()

