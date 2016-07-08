include_directories(Model Solver Utility)

set(DSP_SRC_FILES
	Model/DecData.cpp
	Model/DecDetModel.cpp
	Model/DecTssModel.cpp
	Model/DetModel.cpp
	Model/StoModel.cpp
	Model/TssModel.cpp
	Solver/DecSolver.cpp
	Solver/DspDriver.cpp
	Solver/Benders/BdDriver.cpp
	Solver/Benders/BdMaster.cpp
	Solver/Benders/BdMW.cpp
	Solver/Benders/BdSub.cpp
	Solver/Benders/BdWorker.cpp
	Solver/Benders/SCIPconshdlrBenders.cpp
	Solver/Benders/SCIPconshdlrBendersWorker.cpp
	Solver/Deterministic/DeDriver.cpp
	Solver/DualDecomp/DdDriver.cpp
	Solver/DualDecomp/DdMaster.cpp
	Solver/DualDecomp/DdMasterAtr.cpp
	Solver/DualDecomp/DdMasterDsb.cpp
	Solver/DualDecomp/DdMasterReg.cpp
	Solver/DualDecomp/DdMasterSubgrad.cpp
	Solver/DualDecomp/DdMasterSync.cpp
	Solver/DualDecomp/DdMasterTr.cpp
	Solver/DualDecomp/DdMW.cpp
	Solver/DualDecomp/DdMWAsync.cpp
	Solver/DualDecomp/DdMWPara.cpp
	Solver/DualDecomp/DdMWSync.cpp
	Solver/DualDecomp/DdSub.cpp
	Solver/DualDecomp/DdWorker.cpp
	Solver/DualDecomp/DdWorkerCG.cpp
	Solver/DualDecomp/DdWorkerCGBd.cpp
	Solver/DualDecomp/DdWorkerLB.cpp
	Solver/DualDecomp/DdWorkerUB.cpp
	Solver/DualDecomp/SCIPconshdlrBendersDd.cpp
	SolverInterface/OoqpEps.cpp
	SolverInterface/OoqpStatus.cpp
	SolverInterface/SCIPbranchruleLB.cpp
	SolverInterface/SolverInterfaceClp.cpp
	SolverInterface/SolverInterfaceOoqp.cpp
	SolverInterface/SolverInterfaceOsi.cpp
	SolverInterface/SolverInterfaceScip.cpp
	SolverInterface/SolverInterfaceSpx.cpp
#	SolverInterface/SolverInterfaceCpx.cpp
	Utility/DspMpi.cpp
	Utility/DspParams.cpp)

if (PYDSP_COMPILE)
	set(DSP_SRC_FILES ${DSP_SRC_FILES} PyDsp.cpp)
endif()

add_library(Dsp SHARED
	DspApiEnv.cpp
	DspCInterface.cpp
	${DSP_SRC_FILES})
set(DSP_LIBS
	${SMI_DIR}/lib/libSmi.a
	${SMI_DIR}/lib/libClp.a
	${SMI_DIR}/lib/libOsiClp.a
	${SMI_DIR}/lib/libOsiSpx.a
#	${SMI_DIR}/lib/libOsiCpx.a
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
	${OOQP_DIR}/lib/libooqpbase.a)
if (HAS_MA57)
	set(DSP_LIBS ${DSP_LIBS} ${MA57_DIR}/lib/libma57.a ${METIS_DIR}/libmetis.a)
	message(STATUS "${Blue}DSP will be compiled with MA57 (libma57.a, libmetis.a)${ColourReset}")
else()
	set(DSP_LIBS ${DSP_LIBS} ${MA27_DIR}/lib/libma27.a)
	message(STATUS "${Blue}DSP will be compiled with MA27. Note that DSP can be compiled with MA57 with flag -DHAS_MA57=ON.${ColourReset}")
endif()
set(DSP_LIBS ${DSP_LIBS} 
	gfortran 
	${MPI_CXX_LIBRARIES}
	${BLAS_LIBRARIES}
	${LAPACK_LIBRARIES}
	${ZLIB_LIBRARIES}
	${BZIP2_LIBRARIES})

target_link_libraries(Dsp ${DSP_LIBS})

# INSTALL
install(TARGETS Dsp LIBRARY DESTINATION ${CMAKE_SOURCE_DIR}/lib)

# some message
install(CODE "message(\"\")")
install(CODE "message(\"${Blue}** Please add the following line to the .bash_profile or .bashrc in your home directory:${ColourReset}\")")
if(APPLE)
	install(CODE "message(\"${Red}    export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:${CMAKE_SOURCE_DIR}/lib${ColourReset}\")")
else()
	install(CODE "message(\"${Red}    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${CMAKE_SOURCE_DIR}/lib${ColourReset}\")")
endif()
if(PYDSP_COMPILE)
	if(APPLE)
		install(CODE "execute_process(COMMAND ln -s -f libDsp.dylib lib/_PyDsp.dylib)")
	else()
		install(CODE "execute_process(COMMAND ln -s -f libDsp.so lib/_PyDsp.so)")
	endif()
	install(CODE "execute_process(COMMAND mv src/PyDsp.py lib/PyDsp.py)")
	install(CODE "message(\"${Red}    export PYTHONPATH=$PYTHONPATH:${CMAKE_SOURCE_DIR}/lib${ColourReset}\")")
endif()
install(CODE "message(\"\")")
install(CODE "message(\"** Please install Julia from http://julialang.org/downloads/\")")
install(CODE "message(\"   and add four Julia packages by typing the following in the Julia command-line tool.\")")
install(CODE "message(\"    Pkg.add(\\\"MPI\\\");\")")
install(CODE "message(\"    Pkg.add(\\\"JuMP\\\");\")")
install(CODE "message(\"    Pkg.clone(\\\"https://github.com/kibaekkim/StochJuMP.jl.git\\\");\")")
install(CODE "message(\"    Pkg.clone(\\\"https://github.com/Argonne-National-Laboratory/DSPsolver.jl.git\\\");\")")
install(CODE "message(\"\")")
install(CODE "message(\"** You can test examples:\")")
install(CODE "message(\"    julia examples/farmer/farmer_run.jl BD\")")
install(CODE "message(\"    julia examples/farmer/farmer_run.jl DE\")")
install(CODE "message(\"    mpirun -n 1 julia examples/farmer/farmer_run.jl DD\")")
install(CODE "message(\"\")")
