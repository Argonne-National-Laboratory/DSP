# Julia script of reproducing results for Table 5 in MPC paper.
# Kibaek Kim - ANL MCS 2015

using StochJuMP, DSPsolver

# Algorithmic parameters
nScenarios = int(ARGS[1]); # number of scenarios
prefix     = ARGS[2];

# Unit commitment model
@time include("suc_mod.jl")

DSPsolver.loadProblem(m);

# DSP options
DSPsolver.setLogLevel(3)
DSPsolver.setWallLimit(21600)
DSPsolver.setScipLimitsGap(0.0001)

# SOLVE
DSPsolver.solve(DSP_SOLVER_DE);

# Write results
f = open("output/$prefix.csv", "w");
primal = DSPsolver.getPrimalBound();
dual = DSPsolver.getDualBound();
iter = DSPsolver.getNumIterations();
nodes = DSPsolver.getNumNodes();
println(f, primal, ",", dual, ",", DSPsolver.getSolutionTime(), iter, nodes);
close(f);

println("Solution status: ", DSPsolver.getSolutionStatus());
println("Primal Bound   : ", primal);
println("Dual Bound     : ", dual);
println("Iterations     : ", iter);
println("Nodes          : ", nodes);

