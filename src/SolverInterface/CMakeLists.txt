##################################################
# add coin headers
##################################################
include_directories(
	"${DEPEND_DIR}/include/coin"
	${CPLEX_INC_DIR}
)

set(DSP_SI_FILES)
set(DSP_SI_DEPS Osi Clp)

##################################################
# OOQP interface
##################################################

if(MA27LIB)
	set(DSP_SI_FILES ${DSP_SI_FILES}
		OoqpEps.cpp
		OoqpStatus.cpp
		OsiOoqpSolverInterface.cpp)
	set(DSP_SI_DEPS ${DSP_SI_DEPS} OOQP)
endif()

##################################################
# SCIP interface
##################################################

if(SCIPLIB)
	set(DSP_SI_FILES ${DSP_SI_FILES}
		SCIPbranchruleLB.cpp
		OsiScipSolverInterface.cpp)
endif()

add_library(DspSi OBJECT ${DSP_SI_FILES})
add_dependencies(DspSi ${DSP_SI_DEPS})
set_property(TARGET DspSi PROPERTY POSITION_INDEPENDENT_CODE ON)
