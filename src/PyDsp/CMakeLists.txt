##################################################
# PyDsp library: Python interface
##################################################
add_library(PyDsp SHARED $ENV{PWD}/src/PyDsp/PyDsp.cpp)
target_link_libraries(PyDsp Dsp)
#target_link_libraries(PyDsp Dsp python)
install(TARGETS PyDsp LIBRARY DESTINATION $ENV{PWD}/lib)
