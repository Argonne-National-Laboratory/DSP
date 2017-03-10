##################################################
# REQUIRED PACKAGES
##################################################
find_package(BLAS REQUIRED)

##################################################
# Build and install external packages
##################################################
include(ExternalProject)

find_library(MA27LIB ma27 PATHS ${MA27LIB_DIR})
if(${MA27LIB} MATCHES "MA27LIB-NOTFOUND")
	message(FATAL_ERROR "Library ma27 cannot be found.")
endif()

set(EXTRA_DIR      "${CMAKE_SOURCE_DIR}/extra")
set(COINUTILS_DIR  "${EXTRA_DIR}/CoinUtils")
set(OSI_DIR        "${EXTRA_DIR}/Osi")
set(SMI_DIR        "${EXTRA_DIR}/Smi")
set(ALPS_DIR       "${EXTRA_DIR}/Alps")
set(BCPS_DIR       "${EXTRA_DIR}/Bcps")
set(OOQP_DIR       "${EXTRA_DIR}/OOQP")
set(COIN_INC_DIR   "${DEPEND_DIR}/include/coin")
set(CPLEX_LIB      "-L${CPLEX_LIB_DIR} -lcplex -lm -lpthread")
set(COINUTILS_LIB  "-L${DEPEND_DIR}/lib -lCoinUtils")
set(OSI_LIB        "-L${DEPEND_DIR}/lib -lOsi")
set(ALPS_LIB       "-L${DEPEND_DIR}/lib -lAlps")
set(OOQP_HSL_LIB   "MA27LIB=${MA27LIB}")
set(OOQP_BLAS      "BLAS=${BLAS_LIBRARIES}")
string(REPLACE ";" "\ " OOQP_BLAS "${OOQP_BLAS}")

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
		ADD_CXXFLAGS=-fPIC
		CXX=${CMAKE_CXX_COMPILER}
	BUILD_COMMAND make
	INSTALL_COMMAND make install)

# ExternalProject_Add(Bcps
# 	DEPENDS Alps
# 	SOURCE_DIR ${BCPS_DIR}
# 	CONFIGURE_COMMAND ${BCPS_DIR}/configure
# 		--prefix=${DEPEND_DIR}
# 		--enable-static 
# 		--disable-shared
# 		--with-coinutils-lib=${COINUTILS_LIB}
# 		--with-coinutils-incdir=${COIN_INC_DIR}
# 		--with-osi-lib=${OSI_LIB}
# 		--with-osi-incdir=${COIN_INC_DIR}
# 		--with-alps-lib=${ALPS_LIB}
# 		--with-alps-incdir=${COIN_INC_DIR}
# 		ADD_CXXFLAGS=-fPIC
# 		CXX=${CMAKE_CXX_COMPILER}
# 	BUILD_COMMAND make
# 	INSTALL_COMMAND make install)

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
