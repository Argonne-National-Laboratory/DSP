# This file defines the user-specific settings.

set(MA27LIB_DIR   "")
set(CPLEX_LIB_DIR "")
set(CPLEX_INC_DIR "")
set(SCIP_LIB_DIR  "")
set(SCIP_INC_DIR  "")
set(DEPEND_DIR    $ENV{PWD})

# Please change OFF to ON once the settings are provided.
set(USER_SETTINGS OFF)

if(NOT ${USER_SETTINGS})
	message(FATAL_ERROR "Please complete the user-specific settings in UserConfig.cmake")
endif()
