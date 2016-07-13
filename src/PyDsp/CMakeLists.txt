##################################################
# PyDsp library: Python interface
##################################################
add_library(PyDsp SHARED PyDsp.cpp)
target_link_libraries(PyDsp Dsp)
