# Julia script for DSP (Deterministic Equivalent Form)
# Kibaek Kim - ANL MCS 2014

# DSP solver
include("../julia/Dsp.jl")

# create DSP environment
dsp = DSP()

# read problem from SMPS
readSmps(dsp, "./smps/dcap233_200")

# set print level
setLogLevel(dsp, 1)

# set wallclock limit (sec)
setWallLimit(dsp, 10)

# set node limit
#setNodeLimit(dsp, 10)

# solve deterministic equivalent formulation
solveDe(dsp)

# print results
println("Solution Time   : ", getSolutionTime(dsp))
println("Solution Status : ", getSolutionStatus(dsp))
println("Primal Bound    : ", getPrimalBound(dsp))
println("Dual Bound      : ", getDualBound(dsp))
println("Iterations      : ", getNumIterations(dsp))
println("Number of Nodes : ", getNumNodes(dsp))

