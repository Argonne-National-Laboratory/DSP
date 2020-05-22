if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
	message(STATUS "Extracting library and header information by calling '${CMAKE_Fortran_COMPILER} -v'...")
	execute_process(COMMAND "${CMAKE_Fortran_COMPILER}" "-v" ERROR_VARIABLE
    	GFORTRAN_VERBOSE_STR RESULT_VARIABLE FLAG)

	message(STATUS "'${CMAKE_Fortran_COMPILER} -v' returned:")
	message("${GFORTRAN_VERBOSE_STR}")

	# Detect gfortran version
	string(REGEX MATCH "gcc version [^\t\n ]+" GFORTRAN_VER_STR "${GFORTRAN_VERBOSE_STR}")
	string(REGEX REPLACE "gcc version ([^\t\n ]+)" "\\1" GFORTRAN_VERSION_STRING "${GFORTRAN_VER_STR}")
	message(STATUS "Detected gfortran version ${GFORTRAN_VERSION_STRING}")
	unset(GFORTRAN_VER_STR)
  
	set(MATCH_REGEX "[^\t\n ]+[\t\n ]+")
	set(REPLACE_REGEX "([^\t\n ]+)")
  
	# Find architecture for compiler
	string(REGEX MATCH "Target: [^\t\n ]+"
	  GFORTRAN_ARCH_STR "${GFORTRAN_VERBOSE_STR}")
	message(STATUS "Architecture string: ${GFORTRAN_ARCH_STR}")
	string(REGEX REPLACE "Target: ([^\t\n ]+)" "\\1"
	  GFORTRAN_ARCH "${GFORTRAN_ARCH_STR}")
	message(STATUS "Detected gfortran architecture: ${GFORTRAN_ARCH}")
	unset(GFORTRAN_ARCH_STR)
  
	# Find install prefix, if it exists; if not, use default
	string(REGEX MATCH  "--prefix=[^\t\n ]+[\t\n ]+"
	  GFORTRAN_PREFIX_STR "${GFORTRAN_VERBOSE_STR}")
	if(NOT GFORTRAN_PREFIX_STR)
	  message(STATUS "Detected default gfortran prefix")
	  set(GFORTRAN_PREFIX_DIR "/usr/local") # default prefix for gcc install
	else()
	  string(REGEX REPLACE "--prefix=([^\t\n ]+)" "\\1"
		GFORTRAN_PREFIX_DIR "${GFORTRAN_PREFIX_STR}")
	endif()
	message(STATUS "Detected gfortran prefix: ${GFORTRAN_PREFIX_DIR}")
	unset(GFORTRAN_PREFIX_STR)
  
	# Find install exec-prefix, if it exists; if not, use default
	string(REGEX MATCH "--exec-prefix=[^\t\n ]+[\t\n ]+" "\\1"
	  GFORTRAN_EXEC_PREFIX_STR "${GFORTRAN_VERBOSE_STR}")
	if(NOT GFORTRAN_EXEC_PREFIX_STR)
	  message(STATUS "Detected default gfortran exec-prefix")
	  set(GFORTRAN_EXEC_PREFIX_DIR "${GFORTRAN_PREFIX_DIR}")
	else()
	  string(REGEX REPLACE "--exec-prefix=([^\t\n ]+)" "\\1"
		GFORTRAN_EXEC_PREFIX_DIR "${GFORTRAN_EXEC_PREFIX_STR}")
	endif()
	message(STATUS "Detected gfortran exec-prefix: ${GFORTRAN_EXEC_PREFIX_DIR}")
	UNSET(GFORTRAN_EXEC_PREFIX_STR)

	string(REGEX MATCH "--libdir=[^\t\n ]+"
    GFORTRAN_LIB_DIR_STR "${GFORTRAN_VERBOSE_STR}")
	if(NOT GFORTRAN_LIB_DIR_STR)
		message(STATUS "Found --libdir flag -- not found")
		message(STATUS "Using default gfortran library & include directory paths")
		set(GFORTRAN_LIBRARIES_DIR
			"${GFORTRAN_EXEC_PREFIX_DIR}/lib/gcc/${GFORTRAN_ARCH}/${GFORTRAN_VERSION_STRING}")
  	else(NOT GFORTRAN_LIB_DIR_STR)
    	message(STATUS "Found --libdir flag -- yes")
    	string(REGEX REPLACE "--libdir=([^\t\n ]+)" "\\1"
    		GFORTRAN_LIBRARIES_DIR "${GFORTRAN_LIB_DIR_STR}")
	endif(NOT GFORTRAN_LIB_DIR_STR)
	message(STATUS "gfortran libraries path: ${GFORTRAN_LIBRARIES_DIR}")
	unset(GFORTRAN_LIB_DIR_STR)
endif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")

find_library(GFORTRANLIB REQUIRED NAMES gfortran libgfortran 
	PATHS ENV ${LIBRARY_PATH_VAR_NAME} ${GFORTRAN_LIBRARIES_DIR})