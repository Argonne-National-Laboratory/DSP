##################################################
# Compile subdirectories
##################################################

add_subdirectory(Model)
add_subdirectory(Solver)
add_subdirectory(SolverInterface)
add_subdirectory(TreeSearch)
add_subdirectory(Utility)

##################################################
# Shared library
##################################################

# add_library(DspCore ${DSP_SRC_FILES})
# set_property(TARGET DspCore PROPERTY POSITION_INDEPENDENT_CODE ON)
# add_library(Dsp SHARED DspApiEnv.cpp DspCInterface.cpp)
# add_executable(runDsp dsp.cpp DspApiEnv.cpp DspCInterface.cpp)
# add_dependencies(Dsp Alps OOQP)
# add_dependencies(runDsp Alps OOQP)

# if (APPLE)
# 	target_link_libraries(Dsp DspCore ${DSP_LIBS})
# 	target_link_libraries(runDsp DspCore ${DSP_LIBS})
# else()
# 	target_link_libraries(Dsp DspCore ${DSP_LIBS} "-Wl,--no-undefined")
# 	target_link_libraries(runDsp DspCore ${DSP_LIBS} "-Wl,--no-undefined")
# endif(APPLE)

# ##################################################
# # INSTALL
# ##################################################
# install(TARGETS Dsp LIBRARY DESTINATION $ENV{PWD}/lib)
# install(TARGETS runDsp DESTINATION $ENV{PWD}/bin)

add_library(DspApi OBJECT DspApiEnv.cpp)
add_dependencies(DspApi DspModel DspBaseSolver DspUtil)
