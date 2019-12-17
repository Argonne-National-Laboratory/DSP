set(DSP_BD_FILES
    BdDriver.cpp
    BdDriverSerial.cpp
    BdMaster.cpp
    BdMW.cpp
    BdMWSerial.cpp
    BdSub.cpp
    BdWorker.cpp
    SCIPconshdlrBenders.cpp)

if (MPI_CXX_FOUND)
    set(DSP_BD_FILES ${DSP_BD_FILES}
        BdDriverMpi.cpp
        BdMWMpi.cpp
        SCIPconshdlrBendersWorker.cpp)
endif(MPI_CXX_FOUND)

add_library(DspBd OBJECT ${DSP_BD_FILES})
add_dependencies(DspBd DspDd DspModel DspSi DspUtil)
set_property(TARGET DspBd PROPERTY POSITION_INDEPENDENT_CODE ON)
