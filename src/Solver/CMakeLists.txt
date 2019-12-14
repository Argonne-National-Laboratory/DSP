# Compile deterministic equivalent solver
add_subdirectory(Deterministic)

# Compile Benders solver
if(SCIPLIB)
    add_subdirectory(Benders)
endif(SCIPLIB)

# Compile dual decomposition solver
add_subdirectory(DualDecomp)

# Compile Dantzig-Wolfe decomposition solver
add_subdirectory(DantzigWolfe)

# # Compile Dsp base solver
add_library(DspBaseSolver OBJECT DecSolver.cpp)
add_dependencies(DspBaseSolver DspSi Osi Clp)
