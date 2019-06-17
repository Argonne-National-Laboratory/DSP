# This file defines the user-specific settings.

# Setting CPLEX paths is required.
# The others (ma27 and scip) are optional. For example, OOQP will be
# disabled if MA27LIB_DIR is not set.

set(MA27LIB_DIR   "")
set(CPLEX_LIB_DIR "")
set(CPLEX_INC_DIR "")
set(SCIP_DIR      "")
set(SCIP_LIB_DIR  "")
set(SPX_DIR       "")

# Please change OFF to ON once the settings are provided.
set(USER_SETTINGS OFF)

if(NOT ${USER_SETTINGS})
	message(FATAL_ERROR "Please complete the user-specific settings in UserConfig.cmake")
endif()
set(DEPEND_DIR    $ENV{PWD})
