include_directories(Model Solver Utility)

set(DSP_SRC_FILES
	Solver/OoqpEps.cpp
	Solver/OoqpStatus.cpp
	Solver/SCIPconshdlrBenders.cpp
	Solver/SCIPconshdlrBendersDd.cpp
	Solver/SCIPbranchruleLB.cpp
	Solver/SolverInterfaceOsi.cpp
	Solver/SolverInterfaceClp.cpp
	Solver/SolverInterfaceSpx.cpp
	Solver/SolverInterfaceScip.cpp
	Solver/SolverInterfaceOoqp.cpp
	Solver/TssSolver.cpp
	Solver/TssDe.cpp
	Solver/TssBd.cpp
	Solver/TssBdSub.cpp
	Solver/TssDdMpi.cpp
	Solver/TssDdPrimalMaster.cpp
	Solver/TssDdMasterSubgrad.cpp
	Solver/TssDdSub.cpp
	Solver/TssEval.cpp
	Model/StoModel.cpp
	Model/TssModel.cpp
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
install(CODE "message(\"   and add four Julia pakcages by typing the folloing in the Julia command-line tool.\")")
install(CODE "message(\"    Pkg.add(\\\"MPI\\\");\")")
install(CODE "message(\"    Pkg.add(\\\"JuMP\\\");\")")
install(CODE "message(\"    Pkg.clone(\\\"https://github.com/kibaekkim/StochJuMP.jl.git\\\");\")")
install(CODE "message(\"    Pkg.clone(\\\"https://github.com/kibaekkim/DSPsolver.jl.git\\\");\")")
install(CODE "message(\"\")")
install(CODE "message(\"** You can test examples:\")")
install(CODE "message(\"    julia farmer/farmer_run.jl BD\")")
install(CODE "message(\"    julia farmer/farmer_run.jl DE\")")
install(CODE "message(\"    mpirun -n 1 julia farmer/farmer_run.jl DD\")")
install(CODE "message(\"\")")
