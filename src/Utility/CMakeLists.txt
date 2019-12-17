set(DSP_UTIL_FILES
    DspParams.cpp
)

##################################################
# MPI settings
##################################################

if (MPI_CXX_FOUND)
	set(DSP_UTIL_FILES ${DSP_UTIL_FILES}
		DspMpi.cpp)
endif(MPI_CXX_FOUND)

add_library(DspUtil OBJECT ${DSP_UTIL_FILES})
add_dependencies(DspUtil Osi)
set_property(TARGET DspUtil PROPERTY POSITION_INDEPENDENT_CODE ON)
