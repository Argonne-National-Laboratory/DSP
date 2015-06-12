# Julia script of reproducing results for Table 4 and Figure 5 in MPC paper.
# Kibaek Kim - ANL MCS 2015

using StochJuMP, MPI, DSP

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
DSP.loadProblem(m);

# DSP options
DSP.setLogLevel(1)
DSP.setWallLimit(21600);
DSP.setScipLimitsTime(300);

# DD specific options
DSP.setDdStoppingTolerance(1.0e-4);
DSP.setDdMasterSolver(solver_type);
DSP.setDdMasterNumCutsPerIter(1);
DSP.setDdAddFeasCuts(cuts);
DSP.setDdAddOptCuts(cuts);
DSP.setDdEvalUb(1);

# SOLVE
DSP.solve(DSP_SOLVER_DSP);

if MPI.Comm_rank(MPI.COMM_WORLD) == 0

	# Write results
	f = open("output/$prefix.csv", "w");
	primal = DSP.getPrimalBound();
	dual = DSP.getDualBound();
	println(f, DSP.getNumIterations(), ",", primal, ",", dual, ",", DSP.getSolutionTime());
	close(f);

	if nScenarios == 8
	        writecsv("output/$prefix\_primals.csv", DSP.getDdPrimalBounds());
        	writecsv("output/$prefix\_duals.csv", DSP.getDdDualBounds());
	end

        println("Solution status: ", DSP.getSolutionStatus());
        println("Primal Bound   : ", DSP.getPrimalBound());
        println("Dual Bound     : ", DSP.getDualBound());
end
MPI.Finalize()

