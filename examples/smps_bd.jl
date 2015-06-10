# Julia script for DSP (Benders Decomposition)
# Kibaek Kim - ANL MCS 2014
# Note:
#  1. Integrality restrictions in the second stage will be relaxed.
#  2. This supports only OpenMP (not MPI yet). 

# DSP solver
include("../julia/Dsp.jl")

# create DSP environment
dsp = DSP()

# set print level
setLogLevel(dsp, 1)

# set number of cores
setNumCores(dsp, 1)

# augment scenarios
#augScenarios = [0,1];
#setBdAugScenarios(dsp, 2, augScenarios)

setBendersAggressive(dsp,-200000);

# number of cuts added per iteration
nCuts = 1;

println("\n-----------------")
println("-- dcap233_200 --")
println("-----------------\n")

# read problem from SMPS
readSmps(dsp, "./smps/dcap233_200")

# solve problem using Benders decomposition
solveBd(dsp, nCuts)

# print some results
println("Solution time:   ", getSolutionTime(dsp))
println("Solution status: ", getSolutionStatus(dsp))
println("Objective value: ", getObjValue(dsp))
println("Iterations:      ", getNumIterations(dsp))
println("Number of nodes: ", getNumNodes(dsp))

println("\n-------------------------------")
println("-- Release the current model --")
println("-------------------------------")

# free model to read a new model
freeModel(dsp)

println("\n-----------------")
println("-- dcap233_500 --")
println("-----------------\n")

# read problem from SMPS
readSmps(dsp, "./smps/dcap233_500")

# solve new model
solveBd(dsp, nCuts)

# print some results
println("Solution time:   ", getSolutionTime(dsp))
println("Solution status: ", getSolutionStatus(dsp))
println("Objective value: ", getObjValue(dsp))
println("Iterations:      ", getNumIterations(dsp))
println("Number of nodes: ", getNumNodes(dsp))

println("\n----------------")
println("-- End of run --")
println("----------------")

