add_library(DspDe OBJECT DeDriver.cpp)
add_dependencies(DspDe DspModel DspSi DspUtil)
set_property(TARGET DspDe PROPERTY POSITION_INDEPENDENT_CODE ON)
