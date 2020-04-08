# This file defines the user-specific settings.

set(MA27LIB_DIR     "/usr/local/lib")
#set(CPLEX_LIB_DIR	"")
#set(CPLEX_INC_DIR	"")
set(CPLEX_LIB_DIR   "/Applications/CPLEX_Studio129/opl/lib/x86-64_osx/static_pic")
set(CPLEX_INC_DIR   "/Applications/CPLEX_Studio129/opl/include/ilcplex")
set(GUROBI_LIB_DIR  "/Library/gurobi901/mac64/lib")
set(GUROBI_INC_DIR  "/Library/gurobi901/mac64/include")
set(GUROBI_VERSION  "90")
set(SCIPOPT_INC_DIR "")
set(SCIPOPT_LIB_DIR "")

# Please change OFF to ON once the settings are provided.
set(USER_SETTINGS ON)

if(NOT ${USER_SETTINGS})
	message(FATAL_ERROR "Please complete the user-specific settings in UserConfig.cmake")
endif()
