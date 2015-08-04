# Julia script of reproducing results for Table 4 and Figure 5 in MPC paper.
# Kibaek Kim - ANL MCS 2015

using StochJuMP, MPI, DSPsolver

# Initialize MPI
MPI.Init()

# Algorithmic parameters
nScenarios  = 1; # number of scenarios
solver_type = 2; # solver type
cuts        = -1; # add cuts?
prefix      = "suc";

if length(ARGS) >= 1
	nPeriods  = int(ARGS[1]);
end
if length(ARGS) >= 2
	limit_minupdown = (int(ARGS[2]) != 0)
else
    limit_minupdown = false;
end

# Unit commitment model
include("suc_mod_dec_time.jl")

# Load data to DSP
# DSPsolver.loadProblem(m);

# DSP options
# DSPsolver.setLogLevel(1)
DSPsolver.setLogLevel(3)
# DSPsolver.setWallLimit(21600);
# DSPsolver.setWallLimit(7200);
DSPsolver.setWallLimit(3600);
# DSPsolver.setScipLimitsTime(300);

# DSPsolver.setScipLimitsGap(0.0001)

# DD specific options
DSPsolver.setDdStoppingTolerance(1.0e-5);
# DSPsolver.setDdStoppingTolerance(1.0e-4);
DSPsolver.setDdMasterSolver(solver_type);
DSPsolver.setDdMasterNumCutsPerIter(0);
DSPsolver.setDdAddFeasCuts(cuts);
DSPsolver.setDdAddOptCuts(cuts);
DSPsolver.setDdEvalUb(-1);

DSPsolver.setDdTrustRegionSize(10)
DSPsolver.setDdDisableTrustRegionDecrease(true)

# SOLVE
DSPsolver.solve(DSP_SOLVER_DD);

if MPI.Comm_rank(MPI.COMM_WORLD) == 0

	# Write results
	f = open("$prefix.csv", "w");
	primal = DSPsolver.getPrimalBound();
	dual = DSPsolver.getDualBound();
	println(f, DSPsolver.getNumIterations(), ",", primal, ",", dual, ",", DSPsolver.getSolutionTime());
	close(f);

	if nScenarios == 8
	        writecsv("$prefix\_primals.csv", DSPsolver.getDdPrimalBounds());
        	writecsv("$prefix\_duals.csv", DSPsolver.getDdDualBounds());
	end

        println("Solution status: ", DSPsolver.getSolutionStatus());
        println("Primal Bound   : ", DSPsolver.getPrimalBound());
        println("Dual Bound     : ", DSPsolver.getDualBound());
end
MPI.Finalize()

