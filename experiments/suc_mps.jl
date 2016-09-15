# Julia script of reproducing results for Table 5 in MPC paper.
# Kibaek Kim - ANL MCS 2016

using StochJuMP, DSPsolver

# Algorithmic parameters
nScenarios = parse(Int,ARGS[1]); # number of scenarios
prefix     = ARGS[2];

# Unit commitment model
@time include("suc_mod.jl")

DSPsolver.loadProblem(m);
DSPsolver.writeMps(prefix);
