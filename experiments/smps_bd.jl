# Julia script of reproducing results for Table 2, Figures 1 and 2 in MPC paper.
# Kibaek Kim - ANL MCS 2015

# Arguments from command line
smps_file   = ARGS[1];      # SMPS file name
prefix      = ARGS[2];      # output prefix

using DSPsolver

# read problem from SMPS
DSPsolver.readSmps(smps_file);

# solve problem using Dual Decomposition
DSPsolver.solve(DSP_SOLVER_BD);

# solution results
primal = DSPsolver.getPrimalBound();
dual   = DSPsolver.getDualBound();
gap    = (primal - dual) / abs(primal + 1.0e-10);

# Print out simple results
@printf("Primal Bound : %+e\n", primal);
@printf("Dual Bound   : %+e\n", dual);
@printf("Gap          : %.2f %%\n", gap * 100);

# Write results
f = open("output/$prefix.csv", "w");
println(f, primal, ",", dual, ",", gap);
close(f);

