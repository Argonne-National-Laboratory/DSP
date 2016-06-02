# Julia script
# Kibaek Kim - ANL/MCS 2016

WORK_DIR  = "../examples/smps"
INST_NAME = "sslp_15_45_5";
smps_file = "$WORK_DIR/$INST_NAME";

using MPI, DSPsolver

# Initialize MPI
MPI.Init();

# read problem from SMPS
DSPsolver.readSmps(smps_file);

# set optionts
DSPsolver.setIntParam("LOG_LEVEL", 1);

# solve problem using extensive form
DSPsolver.solve(DSP_SOLVER_BD_MPI);

# solution results
primal = DSPsolver.getPrimalBound();
dual   = DSPsolver.getDualBound();
gap    = (primal - dual) / abs(primal + 1.0e-10);

# Print out simple results
@printf("Primal Bound : %+e\n", primal);
@printf("Dual Bound   : %+e\n", dual);
@printf("Gap          : %.2f %%\n", gap * 100);

# Finalize MPI
MPI.Finalize();
