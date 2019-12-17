set(DSP_DW_FILES
    DwBranchInt.cpp
    DwBranchNonant.cpp
    DwBranchNonant2.cpp
    DwBranchGenDisj.cpp
    DwModel.cpp
    DwMaster.cpp
    DwBundleDual.cpp
    DwBundleDualSmip.cpp
    DwHeuristic.cpp
    DwSolverSerial.cpp
    DwWorker.cpp
)

if(MPI_CXX_FOUND)
    set(DSP_DW_FILES ${DSP_DW_FILES}
        DwSolverMpi.cpp
        DwWorkerMpi.cpp
    )
endif(MPI_CXX_FOUND)

include_directories(.)

add_library(DspDw OBJECT ${DSP_DW_FILES})
add_dependencies(DspDw DspModel DspSi DspUtil Alps)
set_property(TARGET DspDw PROPERTY POSITION_INDEPENDENT_CODE ON)
