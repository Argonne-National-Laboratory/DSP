##################################################
# Build and install external packages
##################################################
include(ExternalProject)

set(EXTRA_DIR      "${CMAKE_SOURCE_DIR}/extra")
set(COINUTILS_DIR  "${EXTRA_DIR}/CoinUtils")
set(OSI_DIR        "${EXTRA_DIR}/Osi")
set(CLP_DIR        "${EXTRA_DIR}/Clp")
set(CGL_DIR        "${EXTRA_DIR}/Cgl")
set(CBC_DIR        "${EXTRA_DIR}/Cbc")
set(SMI_DIR        "${EXTRA_DIR}/Smi")
set(COIN_INC_DIR   "${DEPEND_DIR}/include/coin")
set(COINUTILS_LIB  "-L${DEPEND_DIR}/lib -lCoinUtils")
set(OSI_LIB        "-L${DEPEND_DIR}/lib -lOsi")
set(CLP_LIB        "-L${DEPEND_DIR}/lib -lOsiClp - lClp")
set(CGL_LIB        "-L${DEPEND_DIR}/lib -lCgl")
set(ALPS_DIR       "${EXTRA_DIR}/Alps")
set(OOQP_DIR       "${EXTRA_DIR}/OOQP")
set(ALPS_LIB       "-L${DEPEND_DIR}/lib -lAlps")
set(OOQP_HSL_LIB   "MA27LIB=${MA27LIB}")
set(OOQP_BLAS      "BLAS=${BLAS_LIBRARIES}")
string(REPLACE ";" "\ " OOQP_BLAS "${OOQP_BLAS}")

set(EXTRA_CXX_FLAGS "-fPIC -Wno-sign-conversion -Wno-float-conversion -Wno-shorten-64-to-32 -Wno-unused-variable -Wno-constant-logical-operand -Wno-dynamic-class-memaccess -Wno-unused-private-field")

##################################################
# Compile OOQP with ma27
##################################################

if(MA27LIB)
	ExternalProject_Add(OOQP
		SOURCE_DIR ${OOQP_DIR}
		CONFIGURE_COMMAND ${OOQP_DIR}/configure
			--prefix=${DEPEND_DIR}
			CXXFLAGS=${EXTRA_CXX_FLAGS}
			FFLAGS=-fPIC
			${OOQP_HSL_LIB}
			${OOQP_BLAS}
		BINARY_DIR ${OOQP_DIR}
		INSTALL_DIR ${OOQP_DIR} 
		BUILD_COMMAND make all_libs
		INSTALL_COMMAND make install_headers install_libs)
else()
	message(STATUS "Could not find ma27 (disabled OOQP solver)")
endif()

##################################################
# Compile CoinUtils
##################################################

ExternalProject_Add(CoinUtils
	SOURCE_DIR ${COINUTILS_DIR}
	CONFIGURE_COMMAND ${COINUTILS_DIR}/configure 
		--prefix=${DEPEND_DIR}
		--enable-static 
		--disable-shared
		ADD_CXXFLAGS=${EXTRA_CXX_FLAGS}
	BUILD_COMMAND make -j1
	INSTALL_COMMAND make -j1 install)

##################################################
# Compile OSI
##################################################

if(CPLEXLIB)
	set(CPLEX_LIB "-L${CPLEX_LIB_DIR} -lcplex -lm -lpthread -ldl")
	set(WITH_CPLEX_LIB "--with-cplex-lib=${CPLEX_LIB}")
	set(WITH_CPLEX_INCDIR "--with-cplex-incdir=${CPLEX_INC_DIR}")
else()
	set(WITH_CPLEX_LIB "")
	set(WITH_CPLEX_INCDIR "")
endif()

if(GUROBILIB)
	get_filename_component(LIBGRB ${GUROBILIB} NAME_WE)
	string(SUBSTRING ${LIBGRB} 3 10 GRB)
	set(GUROBI_LIB "-L${GUROBI_LIB_DIR} -l${GRB} -lm -lpthread")
	set(WITH_GUROBI_LIB "--with-gurobi-lib=${GUROBI_LIB}")
	set(WITH_GUROBI_INCDIR "--with-gurobi-incdir=${GUROBI_INC_DIR}")
else()
	set(WITH_GUROBI_LIB "")
	set(WITH_GUROBI_INCDIR "")
endif()

ExternalProject_Add(Osi
	DEPENDS CoinUtils
	SOURCE_DIR ${OSI_DIR}
	CONFIGURE_COMMAND ${OSI_DIR}/configure
		--prefix=${DEPEND_DIR}
		--enable-static 
		--disable-shared
		--with-coinutils-lib=${COINUTILS_LIB}
		--with-coinutils-incdir=${COIN_INC_DIR}
		${WITH_CPLEX_LIB}
		${WITH_CPLEX_INCDIR}
		${WITH_GUROBI_LIB}
		${WITH_GUROBI_INCDIR}
		--disable-cplex-libcheck
		ADD_CXXFLAGS=${EXTRA_CXX_FLAGS}
	BUILD_COMMAND make -j1
	INSTALL_COMMAND make -j1 install)

##################################################
# Compile Clp
##################################################

if(APPLE)
	set(COIN_BLAS "LDFLAGS=-framework Accelerate")
else(APPLE)
	set(COIN_BLAS "--with-blas-lib=${BLAS_LIBRARIES}")
	string(REPLACE ";" "\ " COIN_BLAS "${COIN_BLAS}")
endif(APPLE)

ExternalProject_Add(Clp
	DEPENDS Osi
	SOURCE_DIR ${CLP_DIR}
	CONFIGURE_COMMAND ${CLP_DIR}/configure
		--prefix=${DEPEND_DIR}
		--enable-static 
		--disable-shared
		${COIN_BLAS}
		--with-coinutils-lib=${COINUTILS_LIB}
		--with-coinutils-incdir=${COIN_INC_DIR}
		--with-osi-lib=${OSI_LIB}
		--with-osi-incdir=${COIN_INC_DIR}
		ADD_CXXFLAGS=${EXTRA_CXX_FLAGS}
	BUILD_COMMAND make -j1
	INSTALL_COMMAND make -j1 install)

##################################################
# Compile CoinSMI
##################################################

ExternalProject_Add(Smi
	DEPENDS Osi
	SOURCE_DIR ${SMI_DIR}
	CONFIGURE_COMMAND ${SMI_DIR}/configure
		--prefix=${DEPEND_DIR}
		--enable-static
		--disable-shared
		--with-coinutils-lib=${COINUTILS_LIB}
		--with-coinutils-incdir=${COIN_INC_DIR}
		--with-osi-lib=${OSI_LIB}
		--with-osi-incdir=${COIN_INC_DIR}
		ADD_CXXFLAGS=${EXTRA_CXX_FLAGS}
	BUILD_COMMAND make -j1
	INSTALL_COMMAND make -j1 install)

##################################################
# Compile CoinAlps
##################################################

ExternalProject_Add(Alps
	DEPENDS Osi
	SOURCE_DIR ${ALPS_DIR}
	CONFIGURE_COMMAND ${ALPS_DIR}/configure
		--prefix=${DEPEND_DIR}
		--enable-static 
		--disable-shared
		--with-coinutils-lib=${COINUTILS_LIB}
		--with-coinutils-incdir=${COIN_INC_DIR}
		--with-osi-lib=${OSI_LIB}
		--with-osi-incdir=${COIN_INC_DIR}
		ADD_CXXFLAGS=${EXTRA_CXX_FLAGS}
	BUILD_COMMAND make -j1
	INSTALL_COMMAND make -j1 install)
