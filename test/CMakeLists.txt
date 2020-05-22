add_test(NAME de_farmer COMMAND $ENV{PWD}/src/runDsp --algo de --smps ${CMAKE_SOURCE_DIR}/examples/smps/farmer --param ${CMAKE_SOURCE_DIR}/test/params_cbc.txt --test -108389.9994043)
add_test(NAME dd_farmer COMMAND $ENV{PWD}/src/runDsp --algo dd --smps ${CMAKE_SOURCE_DIR}/examples/smps/farmer --param ${CMAKE_SOURCE_DIR}/test/params_cbc.txt --test -108389.9994043)
add_test(NAME dw_farmer COMMAND ../src/runDsp --algo dw --smps ${CMAKE_SOURCE_DIR}/examples/smps/farmer --param test/params_cbc.txt --test -108389.9994043)

# if(SCIPLIB)
# endif(SCIPLIB)

# if (MPI_CXX_FOUND)
# endif(MPI_CXX_FOUND)