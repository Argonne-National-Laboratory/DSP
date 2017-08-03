##################################################
# Base settings
##################################################

include_directories(
	Model 
	Solver 
	Utility
	${DEPEND_DIR}/include/coin
)

set(DSP_SRC_FILES
	Model/DecData.cpp
	Model/DecDetModel.cpp
	Model/DecTssModel.cpp
	Model/DetModel.cpp
	Model/StoModel.cpp
	Model/TssModel.cpp
	Solver/DecSolver.cpp
	Solver/DspDriver.cpp
	Solver/Deterministic/DeDriver.cpp
	Solver/DualDecomp/DdDriver.cpp
	Solver/DualDecomp/DdDriverSerial.cpp
	Solver/DualDecomp/DdMaster.cpp
	Solver/DualDecomp/DdMasterAtr.cpp
	Solver/DualDecomp/DdMasterSubgrad.cpp
	Solver/DualDecomp/DdMasterTr.cpp
	Solver/DualDecomp/DdMW.cpp
	Solver/DualDecomp/DdMWSerial.cpp
	Solver/DualDecomp/DdSub.cpp
	Solver/DualDecomp/DdWorker.cpp
	Solver/DualDecomp/DdWorkerCG.cpp
	Solver/DualDecomp/DdWorkerLB.cpp
	Solver/DualDecomp/DdWorkerUB.cpp
	SolverInterface/SolverInterfaceClp.cpp
	SolverInterface/SolverInterfaceOsi.cpp
	Utility/DspParams.cpp)

# set DSP libraries
set(DSP_LIBS
	${DEPEND_DIR}/lib/libSmi.a
	${DEPEND_DIR}/lib/libOsiClp.a
	${DEPEND_DIR}/lib/libClp.a
	${DEPEND_DIR}/lib/libOsi.a
	${DEPEND_DIR}/lib/libCoinUtils.a)

##################################################
# MPI settings
##################################################

if(MPI_CXX_FOUND)
	set(DSP_SRC_FILES ${DSP_SRC_FILES}
		Solver/DualDecomp/DdDriverMpi.cpp
		Solver/DualDecomp/DdMWAsync.cpp
		Solver/DualDecomp/DdMWPara.cpp
		Solver/DualDecomp/DdMWSync.cpp
		Utility/DspMpi.cpp)
	set(DSP_LIBS ${DSP_LIBS} ${MPI_CXX_LIBRARIES})
endif(MPI_CXX_FOUND)

##################################################
# OOQP settings
##################################################

if(NOT ${MA27LIB} MATCHES "MA27LIB-NOTFOUND")
	include_directories(
		${DEPEND_DIR}/include/ooqp
	)
	set(DSP_SRC_FILES ${DSP_SRC_FILES}
		Solver/DualDecomp/DdMasterDsb.cpp
		SolverInterface/OoqpEps.cpp
		SolverInterface/OoqpStatus.cpp
		SolverInterface/SolverInterfaceOoqp.cpp)
	set(DSP_LIBS ${DSP_LIBS} 
		${DEPEND_DIR}/lib/libooqpmehrotra.a
		${DEPEND_DIR}/lib/libooqpgondzio.a
		${DEPEND_DIR}/lib/libooqpgendense.a
		${DEPEND_DIR}/lib/libooqpgensparse.a
		${DEPEND_DIR}/lib/libooqpdense.a
		${DEPEND_DIR}/lib/libooqpsparse.a
		${DEPEND_DIR}/lib/libooqpbase.a
		${MA27LIB})
else()
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS} -DNO_OOQP")
endif()

##################################################
# CPLEX settings
##################################################

if(NOT ${CPLEXLIB} MATCHES "CPLEXLIB-NOTFOUND")
	include_directories(
		${CPLEX_INC_DIR}
	)
	set(DSP_SRC_FILES ${DSP_SRC_FILES}
		SolverInterface/SolverInterfaceCpx.cpp)
	set(DSP_LIBS ${DSP_LIBS} 
		${DEPEND_DIR}/lib/libOsiCpx.a
		${CPLEX_LIB_DIR}/libcplex.a)
else()
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS} -DNO_CPX")
endif()

##################################################
# SCIP settings
##################################################

if(NOT ${SCIPLIB} MATCHES "SCIPLIB-NOTFOUND")
	set(DSP_SRC_FILES ${DSP_SRC_FILES}
		Solver/Benders/BdDriver.cpp
		Solver/Benders/BdDriverSerial.cpp
		Solver/Benders/BdMaster.cpp
		Solver/Benders/BdMW.cpp
		Solver/Benders/BdMWSerial.cpp
		Solver/Benders/BdSub.cpp
		Solver/Benders/BdWorker.cpp
		Solver/Benders/SCIPconshdlrBenders.cpp
		Solver/DualDecomp/DdWorkerCGBd.cpp
		Solver/DualDecomp/SCIPconshdlrBendersDd.cpp
		SolverInterface/SCIPbranchruleLB.cpp
		SolverInterface/SolverInterfaceScip.cpp
		SolverInterface/SolverInterfaceSpx.cpp)
	if (MPI_CXX_FOUND)
		set(DSP_SRC_FILES ${DSP_SRC_FILES}
			Solver/Benders/BdDriverMpi.cpp
			Solver/Benders/BdMWMpi.cpp
			Solver/Benders/SCIPconshdlrBendersWorker.cpp)
	endif(MPI_CXX_FOUND)
	set(DSP_LIBS ${DSP_LIBS} 
		${SCIP_DIR}/lib/libscip.a
		${SCIP_DIR}/lib/libobjscip.a
		${SCIP_DIR}/lib/liblpispx.a
		${SCIP_DIR}/lib/libnlpi.cppad.a
		${SOPLEX_DIR}/lib/libsoplex.a)
else()
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS} -DNO_SCIP")
	message(WARNING "SCIP is required for Benders decomposition.")
	message(WARNING "Benders decomposition is disabled with SCIP library.")
endif()

set(DSP_LIBS ${DSP_LIBS}
	gfortran 
	${LAPACK_LIBRARIES}
	${ZLIB_LIBRARIES}
	${BZIP2_LIBRARIES})

add_library(Dsp SHARED
	DspApiEnv.cpp
	DspCInterface.cpp
	${DSP_SRC_FILES})

target_link_libraries(Dsp ${DSP_LIBS} "-Wl,--no-undefined")

##################################################
# INSTALL
##################################################
install(TARGETS Dsp LIBRARY DESTINATION $ENV{PWD}/lib)

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
