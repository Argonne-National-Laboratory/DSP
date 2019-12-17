set(DSP_DD_FILES
    DdDriver.cpp
    DdDriverSerial.cpp
    DdMaster.cpp
    DdMasterAtr.cpp
    DdMasterSubgrad.cpp
    DdMasterTr.cpp
    DdMW.cpp
    DdMWSerial.cpp
    DdSub.cpp
    DdWorkerCG.cpp
    DdWorkerLB.cpp
    DdWorkerUB.cpp
)

if(SCIPLIB)
    set(DSP_DD_FILES ${DSP_DD_FILES} 
        DdWorkerCGBd.cpp
        SCIPconshdlrBendersDd.cpp
    )
endif(SCIPLIB)

if(MPI_CXX_FOUND)
    set(DSP_DD_FILES ${DSP_DD_FILES}
        DdDriverMpi.cpp
        DdMWAsync.cpp
        DdMWAsyncDyn.cpp
        DdMWPara.cpp
        DdMWSync.cpp
    )
endif(MPI_CXX_FOUND)

add_library(DspDd OBJECT ${DSP_DD_FILES})
add_dependencies(DspDd DspModel DspSi DspUtil)
set_property(TARGET DspDd PROPERTY POSITION_INDEPENDENT_CODE ON)
