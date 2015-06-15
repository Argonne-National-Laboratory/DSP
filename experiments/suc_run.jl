# Julia script of reproducing results for Table 4 and Figure 5 in MPC paper.
# Kibaek Kim - ANL MCS 2015

using StochJuMP, MPI, DSPsolver

# Initialize MPI
MPI.Init()

# Algorithmic parameters
nScenarios  = int(ARGS[1]); # number of scenarios
solver_type = int(ARGS[2]); # solver type
cuts        = int(ARGS[3]); # add cuts?
prefix      = ARGS[4];

# Unit commitment model
include("suc_mod.jl")

# Load data to DSP
DSPsolver.loadProblem(m);

# DSP options
DSPsolver.setLogLevel(1)
DSPsolver.setWallLimit(21600);
DSPsolver.setScipLimitsTime(300);

# DD specific options
DSPsolver.setDdStoppingTolerance(1.0e-4);
DSPsolver.setDdMasterSolver(solver_type);
DSPsolver.setDdMasterNumCutsPerIter(1);
DSPsolver.setDdAddFeasCuts(cuts);
DSPsolver.setDdAddOptCuts(cuts);
DSPsolver.setDdEvalUb(1);

# SOLVE
DSPsolver.solve(DSP_SOLVER_DD);

if MPI.Comm_rank(MPI.COMM_WORLD) == 0

	# Write results
	f = open("output/$prefix.csv", "w");
	primal = DSPsolver.getPrimalBound();
	dual = DSPsolver.getDualBound();
	println(f, DSPsolver.getNumIterations(), ",", primal, ",", dual, ",", DSPsolver.getSolutionTime());
	close(f);

	if nScenarios == 8
	        writecsv("output/$prefix\_primals.csv", DSPsolver.getDdPrimalBounds());
        	writecsv("output/$prefix\_duals.csv", DSPsolver.getDdDualBounds());
	end

        println("Solution status: ", DSPsolver.getSolutionStatus());
        println("Primal Bound   : ", DSPsolver.getPrimalBound());
        println("Dual Bound     : ", DSPsolver.getDualBound());
end
MPI.Finalize()

