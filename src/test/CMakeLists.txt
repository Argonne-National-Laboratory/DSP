# Make test executable
set(TEST_SRC 
    tests-main.cpp
    src/tests-DspCInterface.cpp
    src/tests-TssQCModel.cpp
    src/tests-MIQP.cpp
)

include_directories(.)
add_executable(UnitTests ${TEST_SRC} $<TARGET_OBJECTS:DspObj>)
target_link_libraries(UnitTests ${DSP_LIBS})
