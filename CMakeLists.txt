cmake_minimum_required(VERSION 3.5)

# Should use GNU compilers for unix-type OS, including Mac and CygWin
if(UNIX)
	set(CMAKE_C_COMPILER "gcc")
	set(CMAKE_CXX_COMPILER "g++")
	set(CMAKE_Fortran_COMPILER "gfortran")
endif()

##################################################
# Project name
##################################################
project(DSP)
message(STATUS "The system information: ${CMAKE_SYSTEM}.${CMAKE_SYSTEM_PROCESSOR}")

##################################################
# REQUIRED PACKAGES
##################################################
find_package(LAPACK REQUIRED)
find_package(ZLIB REQUIRED)
find_package(BZip2 REQUIRED)
find_package(MPI)

##################################################
# My include files
##################################################
include(UserConfig.cmake)
set(DEPEND_DIR $ENV{PWD})
include(color_scheme.cmake)

##################################################
# Compiler settings
##################################################

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -std=c++0x")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++0x -DNDEBUG")

if(APPLE)
	set(CMAKE_MACOSX_RPATH OFF)
endif()

if (MPI_CXX_FOUND)
	include_directories(${MPI_CXX_INCLUDE_PATH})
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS} -DDSP_HAS_MPI -DOMPI_IGNORE_CXX_SEEK")
endif()

##################################################
# Add include directories
##################################################
include_directories(src)

##################################################
# Add subdirectory for external packages
##################################################
add_subdirectory(extra)
include_directories("${DEPEND_DIR}/include/coin")
set(DSP_LIBS
	${DEPEND_DIR}/lib/libSmi.a 
	${DEPEND_DIR}/lib/libOsi.a 
	${DEPEND_DIR}/lib/libOsiClp.a 
	${DEPEND_DIR}/lib/libClp.a 
	${DEPEND_DIR}/lib/libAlps.a
	${DEPEND_DIR}/lib/libCoinUtils.a 
)

##################################################
# Set flags, paths, libs for external packages
##################################################

if(MPI_CXX_FOUND)
	set(DSP_LIBS ${DSP_LIBS} ${MPI_CXX_LIBRARIES})
endif()

if(MA27LIB)
	include_directories(${DEPEND_DIR}/include/ooqp)
	set(DSP_LIBS ${DSP_LIBS}
		${DEPEND_DIR}/lib/libooqpmehrotra.a
		${DEPEND_DIR}/lib/libooqpgondzio.a
		${DEPEND_DIR}/lib/libooqpgendense.a
		${DEPEND_DIR}/lib/libooqpgensparse.a
		${DEPEND_DIR}/lib/libooqpdense.a
		${DEPEND_DIR}/lib/libooqpsparse.a
		${DEPEND_DIR}/lib/libooqpbase.a
		${MA27LIB}
		${LAPACK_LIBRARIES})
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDSP_HAS_OOQP")
endif(MA27LIB)

if(CPLEXLIB)
	include_directories(${CPLEX_INC_DIR})
    set(DSP_LIBS ${DSP_LIBS}
        ${DEPEND_DIR}/lib/libOsiCpx.a
        ${CPLEX_LIB_DIR}/libcplex.a
        -lm -lpthread -ldl)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDSP_HAS_CPX")
endif(CPLEXLIB)

if(SCIPLIB)
	include_directories(${SCIP_DIR}/src)
    set(DSP_LIBS ${DSP_LIBS}
        ${SCIP_DIR}/lib/libscip.a
        ${SCIP_DIR}/lib/libobjscip.a
        ${SCIP_DIR}/lib/liblpispx.a
        ${SCIP_DIR}/lib/libnlpi.cppad.a
		${SPX_DIR}/lib/libsoplex.a
		-lreadline)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDSP_HAS_SCIP")
else(SCIPLIB)
	message(${Blue} "SCIP is required for Benders decomposition.")
	message(${Blue} "Benders decomposition is disabled with SCIP library.")
endif(SCIPLIB)

set(DSP_LIBS ${DSP_LIBS}
	${ZLIB_LIBRARIES}
	${BZIP2_LIBRARIES}
	-ldl
	gfortran
)

##################################################
# Add subdirectory for the main source code
##################################################
add_subdirectory(src)

set(BASE_OBJS
	$<TARGET_OBJECTS:DspApi>
	$<TARGET_OBJECTS:DspModel>
	$<TARGET_OBJECTS:DspBaseSolver>
	$<TARGET_OBJECTS:DspSi>
	$<TARGET_OBJECTS:DspUtil>)

##################################################
# Testing
##################################################
enable_testing()
add_executable(test_de
	test/test_de.cpp
	${BASE_OBJS}
	$<TARGET_OBJECTS:DspDe>)
add_dependencies(test_de Alps)
target_link_libraries(test_de ${DSP_LIBS})
add_test(NAME de_farmer COMMAND test_de ${CMAKE_SOURCE_DIR}/examples/smps/farmer -108389.9994043)
add_test(NAME de_dcap233_200 COMMAND test_de ${CMAKE_SOURCE_DIR}/examples/smps/dcap233_200 1834.6550481)
add_test(NAME de_sslp_5_25_50 COMMAND test_de ${CMAKE_SOURCE_DIR}/examples/smps/sslp_5_25_50 -121.6)
add_test(NAME de_sslp_15_45_5 COMMAND test_de ${CMAKE_SOURCE_DIR}/examples/smps/sslp_15_45_5 -262.4)

add_executable(test_dd
	test/test_dd.cpp
	${BASE_OBJS}
	$<TARGET_OBJECTS:DspDd>
	$<TARGET_OBJECTS:DspBd>)
add_dependencies(test_dd Alps)
target_link_libraries(test_dd ${DSP_LIBS})
add_test(NAME dd_farmer COMMAND test_dd ${CMAKE_SOURCE_DIR}/examples/smps/farmer -108389.9994043)
add_test(NAME dd_dcap233_200 COMMAND test_dd ${CMAKE_SOURCE_DIR}/examples/smps/dcap233_200 1834.6550481)
add_test(NAME dd_sslp_5_25_50 COMMAND test_dd ${CMAKE_SOURCE_DIR}/examples/smps/sslp_5_25_50 -121.6)
add_test(NAME dd_sslp_15_45_5 COMMAND test_dd ${CMAKE_SOURCE_DIR}/examples/smps/sslp_15_45_5 -262.4)

add_executable(test_bd
test/test_bd.cpp
	${BASE_OBJS}
	$<TARGET_OBJECTS:DspBd>
	$<TARGET_OBJECTS:DspDd>)
target_link_libraries(test_bd ${DSP_LIBS})
add_test(NAME bd_farmer COMMAND test_bd ${CMAKE_SOURCE_DIR}/examples/smps/farmer -108389.9994043)
add_test(NAME bd_dcap233_200 COMMAND test_bd ${CMAKE_SOURCE_DIR}/examples/smps/dcap233_200 882.6151821832)
add_test(NAME bd_sslp_5_25_50 COMMAND test_bd ${CMAKE_SOURCE_DIR}/examples/smps/sslp_5_25_50 -121.6)
add_test(NAME bd_sslp_15_45_5 COMMAND test_bd ${CMAKE_SOURCE_DIR}/examples/smps/sslp_15_45_5 -265.5686127082)

add_executable(test_dw
	test/test_dw.cpp
    ${BASE_OBJS}
    $<TARGET_OBJECTS:DspDw>
    $<TARGET_OBJECTS:DspTs>)
add_dependencies(test_dw Alps Clp)
target_link_libraries(test_dw ${DSP_LIBS})
add_test(NAME dw_farmer COMMAND test_dw ${CMAKE_SOURCE_DIR}/examples/smps/farmer -108389.9994043)
add_test(NAME dw_dcap233_200 COMMAND test_dw ${CMAKE_SOURCE_DIR}/examples/smps/dcap233_200 1834.6550481)
add_test(NAME dw_sslp_5_25_50 COMMAND test_dw ${CMAKE_SOURCE_DIR}/examples/smps/sslp_5_25_50 -121.6)
add_test(NAME dw_sslp_15_45_5 COMMAND test_dw ${CMAKE_SOURCE_DIR}/examples/smps/sslp_15_45_5 -262.4)

# add_library(Dsp SHARED src/DspCInterface.cpp)
# add_executable(runDsp dsp.cpp src/DspCInterface.cpp)

if (MPI_CXX_FOUND)
	add_executable(test_dd_mpi
		test/test_dd_mpi.cpp
		${BASE_OBJS}
		$<TARGET_OBJECTS:DspDd>
		$<TARGET_OBJECTS:DspBd>)
	add_dependencies(test_dd_mpi Alps)
	target_link_libraries(test_dd_mpi ${DSP_LIBS})
	add_test(NAME dd_mpi_farmer COMMAND mpiexec -np 3 ./test_dd_mpi ${CMAKE_SOURCE_DIR}/examples/smps/farmer -108389.9994043)
	add_test(NAME dd_mpi_dcap233_200 COMMAND mpiexec -np 3 ./test_dd_mpi ${CMAKE_SOURCE_DIR}/examples/smps/dcap233_200 1834.6550481)
	add_test(NAME dd_mpi_sslp_5_25_50 COMMAND mpiexec -np 3 ./test_dd_mpi ${CMAKE_SOURCE_DIR}/examples/smps/sslp_5_25_50 -121.6)
	add_test(NAME dd_mpi_sslp_15_45_5 COMMAND mpiexec -np 3 ./test_dd_mpi ${CMAKE_SOURCE_DIR}/examples/smps/sslp_15_45_5 -262.4)

	add_executable(test_bd_mpi
	test/test_bd_mpi.cpp
		${BASE_OBJS}
		$<TARGET_OBJECTS:DspBd>
		$<TARGET_OBJECTS:DspDd>)
	target_link_libraries(test_bd_mpi ${DSP_LIBS})
	add_test(NAME bd_mpi_farmer COMMAND mpiexec -np 3 ./test_bd_mpi ${CMAKE_SOURCE_DIR}/examples/smps/farmer -108389.9994043)
	add_test(NAME bd_mpi_dcap233_200 COMMAND mpiexec -np 3 ./test_bd_mpi ${CMAKE_SOURCE_DIR}/examples/smps/dcap233_200 882.6151821832)
	add_test(NAME bd_mpi_sslp_5_25_50 COMMAND mpiexec -np 3 ./test_bd_mpi ${CMAKE_SOURCE_DIR}/examples/smps/sslp_5_25_50 -121.6)
	add_test(NAME bd_mpi_sslp_15_45_5 COMMAND mpiexec -np 3 ./test_bd_mpi ${CMAKE_SOURCE_DIR}/examples/smps/sslp_15_45_5 -265.5686127082)
endif(MPI_CXX_FOUND)

# if (APPLE)
# 	target_link_libraries(Dsp ${DSP_LIBS})
# 	# target_link_libraries(runDsp ${DSP_LIBS})
# else()
# 	target_link_libraries(Dsp ${DSP_LIBS} "-Wl,--no-undefined")
# 	# target_link_libraries(runDsp ${DSP_LIBS} "-Wl,--no-undefined")
# endif(APPLE)

##################################################
# INSTALL
##################################################
# install(TARGETS Dsp LIBRARY DESTINATION $ENV{PWD}/lib)
# install(TARGETS runDsp DESTINATION $ENV{PWD}/bin)

# ##################################################
# # Messages
# ##################################################
# install(CODE "message(\"\")")
# install(CODE "message(\"${Blue}** Standalone binary file was installed in $ENV{PWD}/bin${ColourReset}\")")
# install(CODE "message(\"\")")
# install(CODE "message(\"${Blue}** Please add the following line to the .bash_profile or .bashrc in your home directory:${ColourReset}\")")
# set(LIBRARY_PATH_VAR_NAME "LD_LIBRARY_PATH")
# if(APPLE)
# 	set(LIBRARY_PATH_VAR_NAME "DYLD_LIBRARY_PATH")
# else()
# endif()
# if(${PIPS_LIB} MATCHES "PIPS_LIB-NOTFOUND")
# 	install(CODE "message(\"${Red}    export ${LIBRARY_PATH_VAR_NAME}=$${LIBRARY_PATH_VAR_NAME}:$ENV{PWD}/lib${ColourReset}\")")
# else()
# 	install(CODE "message(\"${Red}    export ${LIBRARY_PATH_VAR_NAME}=$${LIBRARY_PATH_VAR_NAME}:$ENV{PWD}/lib:${PIPS_LIB_DIR}${ColourReset}\")")
# endif()
