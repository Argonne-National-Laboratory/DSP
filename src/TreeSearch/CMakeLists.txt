##################################################
# add coin headers
##################################################
include_directories(
    "${DEPEND_DIR}/include/coin"
)

set(DSP_TREESEARCH_FILES
    DspModel.cpp
    DspTreeNode.cpp
)

add_library(DspTs OBJECT ${DSP_TREESEARCH_FILES})
add_dependencies(DspTs Alps Osi)
set_property(TARGET DspTs PROPERTY POSITION_INDEPENDENT_CODE ON)
