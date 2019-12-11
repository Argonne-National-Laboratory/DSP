add_library(DspDe OBJECT DeDriver.cpp)
add_dependencies(DspDe DspModel DspSi DspUtil)