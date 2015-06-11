# Julia script for DSP (Deterministic Equivalent Form)
# Kibaek Kim - ANL MCS 2014

using DSP

DSP.readSmps("./smps/dcap233_200"); # read SMPS files
DSP.setLogLevel(1);                 # set print level
DSP.solve(DSP_SOLVER_DE);           # solver with DE

# print results
println("Solution Time   : ", DSP.getSolutionTime());
println("Solution Status : ", DSP.getSolutionStatus());
println("Primal Bound    : ", DSP.getPrimalBound());
println("Dual Bound      : ", DSP.getDualBound());
println("Iterations      : ", DSP.getNumIterations());
println("Number of Nodes : ", DSP.getNumNodes());

