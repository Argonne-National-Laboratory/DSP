##################################################
# add coin headers
##################################################
include_directories(
    "${DEPEND_DIR}/include/coin"
    ../Utility
)

set(DSP_MODEL_FILES
    BlkModel.cpp
    DecBlkModel.cpp
    DecTssModel.cpp
    DetModel.cpp
    StoModel.cpp
    TssModel.cpp
)

add_library(DspModel OBJECT ${DSP_MODEL_FILES})
# set_property(TARGET DspModel PROPERTY POSITION_INDEPENDENT_CODE ON)
add_dependencies(DspModel Smi)