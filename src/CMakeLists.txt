include_directories(Model Solver Utility)

set(DSP_SRC_FILES
	DspDriver.cpp
	Model/BlkModel.cpp
	Model/DecBlkModel.cpp
	Model/DecTssModel.cpp
	Model/DetModel.cpp
	Model/StoModel.cpp
	Model/TssModel.cpp
	Solver/DantzigWolfe/DwMaster.cpp
	Solver/DantzigWolfe/DwMasterTr.cpp
	Solver/DantzigWolfe/DwMasterTrLight.cpp
	Solver/DantzigWolfe/DwSolverSerial.cpp
	Solver/DantzigWolfe/DwWorker.cpp
	SolverInterface/OoqpEps.cpp
	SolverInterface/OoqpStatus.cpp
	SolverInterface/OsiOoqpSolverInterface.cpp
	TreeSearch/DspModel.cpp
	TreeSearch/DspTreeNode.cpp
	Utility/DspParams.cpp)
if (MPI_CXX_FOUND)
	set(DSP_SRC_FILES
		${DSP_SRC_FILES}
		DspDriverMpi.cpp
		Solver/DantzigWolfe/DwSolverMpi.cpp
		Solver/DantzigWolfe/DwWorkerMpi.cpp
		Utility/DspMpi.cpp)
endif()
add_library(Dsp SHARED
	DspApiEnv.cpp
	DspCInterface.cpp
	${DSP_SRC_FILES})

# add HSL library
if (HAS_MA57)
	set(HSL_LIBS ${MA57_DIR}/lib/libma57.a ${METIS_DIR}/libmetis.a)
	message(STATUS "${Blue}DSP will be compiled with MA57 (libma57.a, libmetis.a)${ColourReset}")
else()
	set(HSL_LIBS ${MA27_DIR}/lib/libma27.a)
	message(STATUS "${Blue}DSP will be compiled with MA27. Note that DSP can be compiled with MA57 with flag -DHAS_MA57=ON.${ColourReset}")
endif()

# set DSP libraries
set(DSP_LIBS
	${SMI_DIR}/lib/libSmi.a
	${SMI_DIR}/lib/libClp.a
	${SMI_DIR}/lib/libOsiClp.a
	${CBC_DIR}/lib/libCgl.a
	${CBC_DIR}/lib/libCbc.a
	${CBC_DIR}/lib/libOsiCbc.a
	${SMI_DIR}/lib/libOsiCpx.a
	${SMI_DIR}/lib/libOsi.a
	${SMI_DIR}/lib/libCoinUtils.a
	${ALPS_DIR}/lib/libAlps.a
	${OOQP_DIR}/lib/libooqpmehrotra.a
	${OOQP_DIR}/lib/libooqpgondzio.a
	${OOQP_DIR}/lib/libooqpgendense.a
	${OOQP_DIR}/lib/libooqpgensparse.a
	${OOQP_DIR}/lib/libooqpdense.a
	${OOQP_DIR}/lib/libooqpsparse.a
	${OOQP_DIR}/lib/libooqpbase.a
	${HSL_LIBS}
	${CPLEXLIB_DIR}/libcplex.a
	gfortran
	${BLAS_LIBRARIES}
	${LAPACK_LIBRARIES}
	${ZLIB_LIBRARIES}
	${BZIP2_LIBRARIES})
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
if (PYDSP_COMPILE)
	add_subdirectory(PyDsp)
endif()

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
if(PYDSP_COMPILE)
	if(APPLE)
		install(CODE "execute_process(COMMAND ln -fs $ENV{PWD}/src/PyDsp/libPyDsp.dylib $ENV{PWD}/lib/_PyDsp.dylib)")
	else()
		install(CODE "execute_process(COMMAND ln -fs $ENV{PWD}/src/PyDsp/libPyDsp.so $ENV{PWD}/lib/_PyDsp.so)")
	endif()
	install(CODE "execute_process(COMMAND cp $ENV{PWD}/src/PyDsp/PyDsp.py $ENV{PWD}/lib/PyDsp.py)")
	install(CODE "message(\"${Red}    export PYTHONPATH=$PYTHONPATH:$ENV{PWD}/lib${ColourReset}\")")
endif()
