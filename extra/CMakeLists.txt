##################################################
# Build and install external packages
##################################################
include(ExternalProject)

find_library(MA27LIB  ma27  PATHS ${MA27LIB_DIR})
find_library(CPLEXLIB cplex PATHS ${CPLEX_LIB_DIR})
find_library(SCIPLIB  scip  PATHS ${SCIP_LIB_DIR})

set(EXTRA_DIR      "${CMAKE_SOURCE_DIR}/extra")
set(COINUTILS_DIR  "${EXTRA_DIR}/CoinUtils")
set(OSI_DIR        "${EXTRA_DIR}/Osi")
set(CLP_DIR        "${EXTRA_DIR}/Clp")
set(SMI_DIR        "${EXTRA_DIR}/Smi")
set(COIN_INC_DIR   "${DEPEND_DIR}/include/coin")
set(COINUTILS_LIB  "-L${DEPEND_DIR}/lib -lCoinUtils")
set(OSI_LIB        "-L${DEPEND_DIR}/lib -lOsi")

##################################################
# Compile OOQP with ma27
##################################################

if(${MA27LIB} MATCHES "MA27LIB-NOTFOUND")
	message(WARNING "Library ma27 cannot be found.")
	message(WARNING "OOQP is disabled.")
else()
	set(OOQP_DIR     "${EXTRA_DIR}/OOQP")
	set(OOQP_HSL_LIB "MA27LIB=${MA27LIB}")
	set(OOQP_BLAS    "BLAS=${BLAS_LIBRARIES}")

	ExternalProject_Add(OOQP
		SOURCE_DIR ${OOQP_DIR}
		CONFIGURE_COMMAND ${OOQP_DIR}/configure
			--prefix=${DEPEND_DIR}
			CXXFLAGS=-fPIC 
			FFLAGS=-fPIC
			${OOQP_HSL_LIB}
			${OOQP_BLAS}
			CXX=${CMAKE_CXX_COMPILER}
		BINARY_DIR ${OOQP_DIR}
		INSTALL_DIR ${OOQP_DIR} 
		BUILD_COMMAND make all_libs
		INSTALL_COMMAND make install_headers install_libs)
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
		ADD_CXXFLAGS=-fPIC
		CXX=${CMAKE_CXX_COMPILER}
	BUILD_COMMAND make
	INSTALL_COMMAND make install)

##################################################
# Compile OSI
##################################################

if(${CPLEXLIB} MATCHES "CPLEXLIB-NOTFOUND")
if(${SPX_DIR} EQUAL "")
	ExternalProject_Add(Osi
		DEPENDS CoinUtils
		SOURCE_DIR ${OSI_DIR}
		CONFIGURE_COMMAND ${OSI_DIR}/configure
			--prefix=${DEPEND_DIR}
			--enable-static 
			--disable-shared
			--with-coinutils-lib=${COINUTILS_LIB}
			--with-coinutils-incdir=${COIN_INC_DIR}
			ADD_CXXFLAGS=-fPIC
			CXX=${CMAKE_CXX_COMPILER}
		BUILD_COMMAND make
		INSTALL_COMMAND make install)
else()
	set(SPX_LIB "-L${SPX_DIR}/lib -lsoplex")
	ExternalProject_Add(Osi
		DEPENDS CoinUtils
		SOURCE_DIR ${OSI_DIR}
		CONFIGURE_COMMAND ${OSI_DIR}/configure
			--prefix=${DEPEND_DIR}
			--enable-static 
			--disable-shared
			--with-coinutils-lib=${COINUTILS_LIB}
			--with-coinutils-incdir=${COIN_INC_DIR}
			--with-soplex-incdir=${SPX_DIR}/src
			--with-soplex-lib=${SPX_LIB}
			ADD_CXXFLAGS=-fPIC
			CXX=${CMAKE_CXX_COMPILER}
		BUILD_COMMAND make
		INSTALL_COMMAND make install)
endif()
else()
	set(CPLEX_LIB "-L${CPLEX_LIB_DIR} -lcplex -lm -lpthread")
	ExternalProject_Add(Osi
		DEPENDS CoinUtils
		SOURCE_DIR ${OSI_DIR}
		CONFIGURE_COMMAND ${OSI_DIR}/configure
			--prefix=${DEPEND_DIR}
			--enable-static 
			--disable-shared
			--with-coinutils-lib=${COINUTILS_LIB}
			--with-coinutils-incdir=${COIN_INC_DIR}
			--with-cplex-lib=${CPLEX_LIB}
			--with-cplex-incdir=${CPLEX_INC_DIR}
			--disable-cplex-libcheck
			ADD_CXXFLAGS=-fPIC
			CXX=${CMAKE_CXX_COMPILER}
		BUILD_COMMAND make
		INSTALL_COMMAND make install)
endif()

##################################################
# Compile Clp
##################################################

if(APPLE)
	set(CLP_BLAS "LDFLAGS=-framework Accelerate")
else(APPLE)
	set(CLP_BLAS "--with-blas-lib=${BLAS_LIBRARIES}")
endif(APPLE)

ExternalProject_Add(Clp
	DEPENDS Osi
	SOURCE_DIR ${CLP_DIR}
	CONFIGURE_COMMAND ${CLP_DIR}/configure
		--prefix=${DEPEND_DIR}
		--enable-static 
		--disable-shared
		${CLP_BLAS}
		--with-coinutils-lib=${COINUTILS_LIB}
		--with-coinutils-incdir=${COIN_INC_DIR}
		--with-osi-lib=${OSI_LIB}
		--with-osi-incdir=${COIN_INC_DIR}
		ADD_CXXFLAGS=-fPIC
		CXX=${CMAKE_CXX_COMPILER}
	BUILD_COMMAND make
	INSTALL_COMMAND make install)

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
		ADD_CXXFLAGS=-fPIC
		CXX=${CMAKE_CXX_COMPILER}
	BUILD_COMMAND make
	INSTALL_COMMAND make install)

