include_directories(Model Solver Utility)

set(DSP_SRC_FILES
	SolverInterface/OoqpEps.cpp
	SolverInterface/OoqpStatus.cpp
	SolverInterface/SCIPconshdlrBenders.cpp
	SolverInterface/SCIPconshdlrBendersDd.cpp
	SolverInterface/SCIPbranchruleLB.cpp
	SolverInterface/SolverInterfaceOsi.cpp
	SolverInterface/SolverInterfaceClp.cpp
	SolverInterface/SolverInterfaceSpx.cpp
	SolverInterface/SolverInterfaceScip.cpp
	SolverInterface/SolverInterfaceOoqp.cpp
	Solver/TssSolver.cpp
	Solver/TssBd.cpp
	Solver/TssBdSub.cpp
	Solver/TssEval.cpp
	Solver/DecSolver.cpp
	Solver/DecDdMpi.cpp
	Solver/DecDdPrimalMaster.cpp
	Solver/DecDdMasterSubgrad.cpp
	Solver/DecDdSub.cpp
	Solver/DecDe.cpp
	Model/StoModel.cpp
	Model/TssModel.cpp
	Model/DecData.cpp
	Model/DecTssModel.cpp
	Model/DecDetModel.cpp
	Model/DetModel.cpp
	Utility/StoMessage.cpp
	Utility/StoUtility.cpp)

add_library(Dsp SHARED
	StoApiEnv.cpp
	StoCInterface.cpp
	${DSP_SRC_FILES})

set(DSP_LIBS
	${SMI_DIR}/lib/libSmi.a
	${SMI_DIR}/lib/libClp.a
	${SMI_DIR}/lib/libOsiClp.a
	${SMI_DIR}/lib/libOsiSpx.a
	${SMI_DIR}/lib/libOsi.a
	${SMI_DIR}/lib/libCoinUtils.a
	${SCIP_DIR}/lib/libscip.a
	${SCIP_DIR}/lib/libobjscip.a
	${SCIP_DIR}/lib/liblpispx.a
	${SCIP_DIR}/lib/libnlpi.cppad.a
	${SOPLEX_DIR}/lib/libsoplex.a
	${OOQP_DIR}/lib/libooqpmehrotra.a
	${OOQP_DIR}/lib/libooqpgondzio.a
	${OOQP_DIR}/lib/libooqpgensparse.a
	${OOQP_DIR}/lib/libooqpsparse.a
	${OOQP_DIR}/lib/libooqpbase.a
	${MA27_DIR}/lib/libma27.a
	${MPI_CXX_LIBRARIES}
	${BLAS_LIBRARIES}
	${LAPACK_LIBRARIES}
	${ZLIB_LIBRARIES}
	${BZIP2_LIBRARIES})

if(APPLE)
	set(DSP_LIBS ${DSP_LIBS} gfortran)
endif()

target_link_libraries(Dsp ${DSP_LIBS})

# INSTALL
install(TARGETS Dsp LIBRARY DESTINATION ${CMAKE_SOURCE_DIR}/lib)

# some message
install(CODE "message(\"\")")
install(CODE "message(\"** Please add the following line to the .bash_profile or .bashrc in your home directory:\")")
if(APPLE)
	install(CODE "message(\"    export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:${CMAKE_SOURCE_DIR}/lib\")")
else()
	install(CODE "message(\"    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${CMAKE_SOURCE_DIR}/lib\")")
endif()
install(CODE "message(\"\")")
install(CODE "message(\"** Please install Julia from http://julialang.org/downloads/\")")
install(CODE "message(\"   and add four Julia packages by typing the following in the Julia command-line tool.\")")
install(CODE "message(\"    Pkg.add(\\\"MPI\\\");\")")
install(CODE "message(\"    Pkg.add(\\\"JuMP\\\");\")")
install(CODE "message(\"    Pkg.clone(\\\"https://github.com/kibaekkim/StochJuMP.jl.git\\\");\")")
install(CODE "message(\"    Pkg.clone(\\\"https://github.com/kibaekkim/DSPsolver.jl.git\\\");\")")
install(CODE "message(\"\")")
install(CODE "message(\"** You can test examples:\")")
install(CODE "message(\"    julia examples/farmer/farmer_run.jl BD\")")
install(CODE "message(\"    julia examples/farmer/farmer_run.jl DE\")")
install(CODE "message(\"    mpirun -n 1 julia examples/farmer/farmer_run.jl DD\")")
install(CODE "message(\"\")")
