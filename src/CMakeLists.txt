# Source codes are stored in variables:
# DSP_SRC

# Include directories are stored in
# DSP_INC_DIR

# Libraries are stored in
# DSP_LIBS

set(DSP_SRC
    DspApiEnv.cpp
    DspCInterface.cpp
    Model/BlkModel.cpp
    Model/DecBlkModel.cpp
    Model/DecTssModel.cpp
    Model/DetModel.cpp
    Model/StoModel.cpp
    Model/TssModel.cpp
    Solver/DecSolver.cpp
    Solver/DantzigWolfe/DwBranchInt.cpp
    Solver/DantzigWolfe/DwBranchNonant.cpp
    Solver/DantzigWolfe/DwBranchNonant2.cpp
    Solver/DantzigWolfe/DwBranchGenDisj.cpp
    Solver/DantzigWolfe/DwModel.cpp
    Solver/DantzigWolfe/DwModelSmip.cpp
    Solver/DantzigWolfe/DwMaster.cpp
    Solver/DantzigWolfe/DwBundleDual.cpp
    Solver/DantzigWolfe/DwBundleDualSmip.cpp
    Solver/DantzigWolfe/DwHeuristic.cpp
    Solver/DantzigWolfe/DwSolverSerial.cpp
    Solver/DantzigWolfe/DwWorker.cpp
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
    Solver/DualDecomp/DdWorkerCG.cpp
    Solver/DualDecomp/DdWorkerLB.cpp
    Solver/DualDecomp/DdWorkerUB.cpp
    TreeSearch/DspModel.cpp
    TreeSearch/DspTreeNode.cpp
    Utility/DspParams.cpp
)
set(DSP_INC_DIR 
    .
    ${DEPEND_DIR}/include/coin)
set(DSP_LIBS 
    ${MPI_CXX_LIBRARIES}
	${DEPEND_DIR}/lib/libSmi.a 
	${DEPEND_DIR}/lib/libOsi.a 
	${DEPEND_DIR}/lib/libOsiClp.a 
	${DEPEND_DIR}/lib/libClp.a 
	${DEPEND_DIR}/lib/libAlps.a
	${DEPEND_DIR}/lib/libCoinUtils.a 
	${ZLIB_LIBRARIES}
	${BZIP2_LIBRARIES}
    ${GFORTRANLIB}
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
    -ldl
)

if(CPLEXLIB)
    set(DSP_SRC 
        ${DSP_SRC}
        Solver/DantzigWolfe/DwBranchInt.cpp
        Solver/DantzigWolfe/DwBranchNonant.cpp
        Solver/DantzigWolfe/DwBranchNonant2.cpp
        Solver/DantzigWolfe/DwBranchGenDisj.cpp
        Solver/DantzigWolfe/DwModel.cpp
        Solver/DantzigWolfe/DwMaster.cpp
        Solver/DantzigWolfe/DwBundleDual.cpp
        Solver/DantzigWolfe/DwBundleDualSmip.cpp
        Solver/DantzigWolfe/DwHeuristic.cpp
        Solver/DantzigWolfe/DwSolverSerial.cpp
        Solver/DantzigWolfe/DwWorker.cpp
    )
    set(DSP_INC_DIR 
        ${DSP_INC_DIR} 
        ${CPLEX_INC_DIR})
    set(DSP_LIBS 
        ${DSP_LIBS}
        ${DEPEND_DIR}/lib/libOsiCpx.a
        ${CPLEX_LIB_DIR}/libcplex.a
        -lm -lpthread -ldl)
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDSP_HAS_CPX")
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DDSP_HAS_CPX")
endif(CPLEXLIB)

if(SCIPLIB)
    set(DSP_SRC 
        ${DSP_SRC}
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
		SolverInterface/OsiScipSolverInterface.cpp
    )
    set(DSP_INC_DIR 
        ${DSP_INC_DIR} 
        ${SCIPOPT_INC_DIR})
    set(DSP_LIBS 
        ${DSP_LIBS}
		${SCIPLIB} -lreadline)
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDSP_HAS_SCIP")
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DDSP_HAS_SCIP")
else(SCIPLIB)
    message(STATUS ${Blue} "SCIP is required for Benders decomposition. Benders decomposition is disabled with SCIP library.${ColourReset}")
endif(SCIPLIB)

if(MA27LIB)
    set(DSP_SRC 
        ${DSP_SRC}
        SolverInterface/OoqpEps.cpp
        SolverInterface/OoqpStatus.cpp
        SolverInterface/OsiOoqpSolverInterface.cpp
    )
    set(DSP_INC_DIR 
        ${DSP_INC_DIR} 
        ${DEPEND_DIR}/include/ooqp)
    set(DSP_LIBS 
        ${DSP_LIBS}
		${DEPEND_DIR}/lib/libooqpmehrotra.a
		${DEPEND_DIR}/lib/libooqpgondzio.a
		${DEPEND_DIR}/lib/libooqpgendense.a
		${DEPEND_DIR}/lib/libooqpgensparse.a
		${DEPEND_DIR}/lib/libooqpdense.a
		${DEPEND_DIR}/lib/libooqpsparse.a
		${DEPEND_DIR}/lib/libooqpbase.a
        ${MA27LIB}
    )
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDSP_HAS_OOQP")
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DDSP_HAS_OOQP")
endif(MA27LIB)

if (MPI_CXX_FOUND)
    set(DSP_SRC 
        ${DSP_SRC}
        Solver/DantzigWolfe/DwSolverMpi.cpp
        Solver/DantzigWolfe/DwWorkerMpi.cpp
        Solver/DualDecomp/DdDriverMpi.cpp
        Solver/DualDecomp/DdMWAsync.cpp
        Solver/DualDecomp/DdMWAsyncDyn.cpp
        Solver/DualDecomp/DdMWPara.cpp
        Solver/DualDecomp/DdMWSync.cpp
        Utility/DspMpi.cpp
    )
    if(SCIPLIB)
        set(DSP_SRC
            ${DSP_SRC}
            Solver/Benders/BdDriverMpi.cpp
            Solver/Benders/BdMWMpi.cpp
            Solver/Benders/SCIPconshdlrBendersWorker.cpp
        )
    endif(SCIPLIB)
    set(DSP_INC_DIR 
        ${DSP_INC_DIR} 
        ${MPI_CXX_INCLUDE_PATH})
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${MPI_CXX_COMPILE_FLAGS} -DDSP_HAS_MPI -DOMPI_IGNORE_CXX_SEEK")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${MPI_CXX_COMPILE_FLAGS} -DDSP_HAS_MPI -DOMPI_IGNORE_CXX_SEEK")
endif(MPI_CXX_FOUND)

include_directories(${DSP_INC_DIR})
add_library(DspObj OBJECT ${DSP_SRC})
set_property(TARGET DspObj PROPERTY POSITION_INDEPENDENT_CODE ON)
add_library(Dsp SHARED $<TARGET_OBJECTS:DspObj>)
add_executable(runDsp dsp.cpp $<TARGET_OBJECTS:DspObj>)

if (APPLE)
    set(CMAKE_EXE_LINKER_FLAGS "-Wl,-undefined,error")
    set(CMAKE_SHARED_LINKER_FLAGS "-Wl,-undefined,error")
else(APPLE)
    set(CMAKE_EXE_LINKER_FLAGS "-Wl,--no-undefined")
    set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined")
endif(APPLE)
target_link_libraries(Dsp ${DSP_LIBS})
target_link_libraries(runDsp ${DSP_LIBS})

##################################################
# INSTALL
##################################################
install(TARGETS Dsp LIBRARY DESTINATION $ENV{PWD}/lib)
install(TARGETS runDsp DESTINATION $ENV{PWD}/bin)

install(CODE "message(\"\")")
install(CODE "message(\"${Blue}** Standalone binary file was installed in $ENV{PWD}/bin${ColourReset}\")")
install(CODE "message(\"\")")
install(CODE "message(\"${Blue}** Please add the following line to the .bash_profile or .bashrc in your home directory:${ColourReset}\")")

##################################################
# Test
##################################################
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../test)
add_executable(test_de ../test/test_de.cpp $<TARGET_OBJECTS:DspObj>)
add_executable(test_dd ../test/test_dd.cpp $<TARGET_OBJECTS:DspObj>)
target_link_libraries(test_de ${DSP_LIBS})
target_link_libraries(test_dd ${DSP_LIBS})
if(SCIPLIB)
    add_executable(test_bd ../test/test_bd.cpp $<TARGET_OBJECTS:DspObj>)
    target_link_libraries(test_bd ${DSP_LIBS})
endif(SCIPLIB)
if(CPLEXLIB)
    add_executable(test_dw ../test/test_dw.cpp $<TARGET_OBJECTS:DspObj>)
    target_link_libraries(test_dw ${DSP_LIBS})
endif(CPLEXLIB)
if(MPI_CXX_FOUND)
    add_executable(test_dd_mpi ../test/test_dd_mpi.cpp $<TARGET_OBJECTS:DspObj>)
    target_link_libraries(test_dd_mpi ${DSP_LIBS})
    if(SCIPLIB)
        add_executable(test_bd_mpi ../test/test_bd_mpi.cpp $<TARGET_OBJECTS:DspObj>)
        target_link_libraries(test_bd_mpi ${DSP_LIBS})
    endif(SCIPLIB)
    if(CPLEXLIB)
        add_executable(test_dw_mpi ../test/test_dw_mpi.cpp $<TARGET_OBJECTS:DspObj>)
        target_link_libraries(test_dw_mpi ${DSP_LIBS})
    endif(CPLEXLIB)
endif(MPI_CXX_FOUND)
