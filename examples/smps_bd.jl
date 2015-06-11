# Julia script for DSP (Benders Decomposition)
# Kibaek Kim - ANL MCS 2014
# Note: Integrality restrictions in the second stage will be relaxed.

using DSP

DSP.readSmps("./smps/dcap233_200"); # read SMPS files
DSP.setLogLevel(1);                 # set print level
DSP.solve(DSP_SOLVER_BD);           # solve with BD

# print some results
println("Solution Time   : ", DSP.getSolutionTime());
println("Solution Status : ", DSP.getSolutionStatus());
println("Primal Bound    : ", DSP.getPrimalBound());
println("Dual Bound      : ", DSP.getDualBound());
println("Iterations      : ", DSP.getNumIterations());

