# Julia script of reproducing results for Table 2, Figures 1 and 2 in MPC paper.
# Kibaek Kim - ANL MCS 2015

smps_file = "./smps/sslp_5_25_50";
prefix    = "sslp_5_25_50";

# Arguments from command line
if length(ARGS) == 2
	smps_file = ARGS[1]; # SMPS file name
	prefix    = ARGS[2]; # output prefix
end

using DSPsolver

# read problem from SMPS
DSPsolver.readSmps(smps_file);

# set optionts
DSPsolver.setLogLevel(1);
DSPsolver.setIntRelax(DSP_SECOND_STAGE);

# solve problem using extensive form
DSPsolver.solve(DSP_SOLVER_DE);

# solution results
primal = DSPsolver.getPrimalBound();
dual   = DSPsolver.getDualBound();
gap    = (primal - dual) / abs(primal + 1.0e-10);

# Print out simple results
@printf("Primal Bound : %+e\n", primal);
@printf("Dual Bound   : %+e\n", dual);
@printf("Gap          : %.2f %%\n", gap * 100);

# Write results
f = open("$prefix.csv", "w");
println(f, primal, ",", dual, ",", gap);
close(f);
