##################################################
# Build and install external packages
##################################################
include(ExternalProject)

set(EXTRA_DIR      "${CMAKE_SOURCE_DIR}/extra")
set(COINUTILS_DIR  "${EXTRA_DIR}/CoinUtils")
set(OSI_DIR        "${EXTRA_DIR}/Osi")
set(CLP_DIR        "${EXTRA_DIR}/Clp")
set(SMI_DIR        "${EXTRA_DIR}/Smi")
set(COIN_INC_DIR   "${DEPEND_DIR}/include/coin-or")
set(COINUTILS_LIB  "-L${DEPEND_DIR}/lib -lCoinUtils")
set(OSI_LIB        "-L${DEPEND_DIR}/lib -lOsi")
set(ALPS_DIR       "${EXTRA_DIR}/Alps")
set(BCPS_DIR       "${EXTRA_DIR}/Bcps")
set(OOQP_DIR       "${EXTRA_DIR}/OOQP")
set(ALPS_LIB       "-L${DEPEND_DIR}/lib -lAlps")
set(OOQP_HSL_LIB   "MA27LIB=${MA27LIB}")
set(OOQP_BLAS      "BLAS=${BLAS_LIBRARIES}")
string(REPLACE ";" "\ " OOQP_BLAS "${OOQP_BLAS}")

##################################################
# Compile OOQP with ma27
##################################################

if(MA27LIB)
	ExternalProject_Add(OOQP
		SOURCE_DIR ${OOQP_DIR}
		CONFIGURE_COMMAND ${OOQP_DIR}/configure
			--prefix=${DEPEND_DIR}
			CXXFLAGS=-fPIC 
			FFLAGS=-fPIC
			${OOQP_HSL_LIB}
			${OOQP_BLAS}
			# CXX=${CMAKE_CXX_COMPILER}
		BINARY_DIR ${OOQP_DIR}
		INSTALL_DIR ${OOQP_DIR} 
		BUILD_COMMAND make all_libs
		INSTALL_COMMAND make install_headers install_libs)
else()
	message(WARNING "Library ma27 cannot be found.")
	message(WARNING "OOQP is disabled.")
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
		# CXX=${CMAKE_CXX_COMPILER}
	BUILD_COMMAND make
	INSTALL_COMMAND make install)

##################################################
# Compile OSI
##################################################

ExternalProject_Add(Glpk
	DOWNLOAD_DIR ${EXTRA_DIR}
	SOURCE_DIR ${EXTRA_DIR}/glpk-4.65
	INSTALL_DIR ${EXTRA_DIR}/glpk-4.65/build
	DOWNLOAD_COMMAND wget http://ftp.gnu.org/gnu/glpk/glpk-4.65.tar.gz &&
		tar -zxvf glpk-4.65.tar.gz && rm glpk-4.65.tar.gz
	CONFIGURE_COMMAND ${EXTRA_DIR}/glpk-4.65/configure --prefix=<INSTALL_DIR> CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} CFLAGS=${EXTRA_CMAKE_INCLUDE_FLAGS} CPPFLAGS=${EXTRA_CMAKE_INCLUDE_FLAGS} LDFLAGS=${EXTRA_CMAKE_LIBRARY_FLAGS}
	BUILD_COMMAND make
	INSTALL_COMMAND make install
)
set(GLPK_LIB "-L${EXTRA_DIR}/glpk-4.65/build/lib -lglpk")
set(GLPK_INCDIR "${EXTRA_DIR}/glpk-4.65/build/include")

if(CPLEXLIB)
	set(CPLEX_LIB "-L${CPLEX_LIB_DIR} -lcplex -lm -lpthread")
	set(WITH_CPLEX_LIB "--with-cplex --with-cplex-lflags=${CPLEX_LIB}")
	set(WITH_CPLEX_INCDIR "--with-cplex-cflags=${CPLEX_INC_DIR}")
else()
	set(WITH_CPLEX_LIB "")
	set(WITH_CPLEX_INCDIR "")
endif()

ExternalProject_Add(Osi
	DEPENDS CoinUtils Glpk
	SOURCE_DIR ${OSI_DIR}
	CONFIGURE_COMMAND ${OSI_DIR}/configure
		--prefix=${DEPEND_DIR}
		--enable-static 
		--disable-shared
		--with-coinutils
		--with-coinutils-lflags=${COINUTILS_LIB}
		--with-coinutils-cflags=-I${COIN_INC_DIR}
		--with-glpk
		--with-glpk-lflags=${GLPK_LIB}
		--with-glpk-cflags=-I${GLPK_INCDIR}
		${WITH_CPLEX_LIB}
		${WITH_CPLEX_INCDIR}
		--disable-cplex-libcheck
		ADD_CXXFLAGS=-fPIC
		# CXX=${CMAKE_CXX_COMPILER}
	BUILD_COMMAND make
	INSTALL_COMMAND make install)

##################################################
# Compile Clp
##################################################

if(APPLE)
	set(CLP_BLAS "LDFLAGS=-framework Accelerate")
else(APPLE)
	set(CLP_BLAS "--with-blas-lib=${BLAS_LIBRARIES}")
    string(REPLACE ";" "\ " CLP_BLAS "${CLP_BLAS}")
endif(APPLE)

ExternalProject_Add(Clp
	DEPENDS Osi
	SOURCE_DIR ${CLP_DIR}
	CONFIGURE_COMMAND ${CLP_DIR}/configure
		--prefix=${DEPEND_DIR}
		--enable-static 
		--disable-shared
		${CLP_BLAS}
		--with-coinutils
		--with-coinutils-lflags=${COINUTILS_LIB}
		--with-coinutils-cflags=-I${COIN_INC_DIR}
		--with-osi
		--with-osi-lflags=${OSI_LIB}
		--with-osi-cflags=-I${COIN_INC_DIR}
		ADD_CXXFLAGS=-fPIC
		# CXX=${CMAKE_CXX_COMPILER}
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
		# --enable-debug
		--disable-shared
		--with-coinutils-lib=${COINUTILS_LIB}
		--with-coinutils-incdir=${COIN_INC_DIR}
		--with-osi-lib=${OSI_LIB}
		--with-osi-incdir=${COIN_INC_DIR}
		ADD_CXXFLAGS=-fPIC
		# CXX=${CMAKE_CXX_COMPILER}
	BUILD_COMMAND make
	INSTALL_COMMAND make install)

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
		ADD_CXXFLAGS=-fPIC
		# CXX=${CMAKE_CXX_COMPILER}
	BUILD_COMMAND make
	INSTALL_COMMAND make install)
