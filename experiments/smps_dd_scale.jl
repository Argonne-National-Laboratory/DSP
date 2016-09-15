# Julia script of reproducing scaling results
# Kibaek Kim - ANL MCS 2016

# Arguments from command line
smps_file   = ARGS[1];            # SMPS file name
wtime       = parse(Int,ARGS[2]); # wall time limit
cut         = parse(Int,ARGS[3]); # use cut?
maxeval     = parse(Int,ARGS[4]); # maximum number of solutions to evaluate
prefix      = ARGS[5];            # output prefix

using MPI, DSPsolver

# initialize MPI
MPI.Init()

# read problem from SMPS
DSPsolver.readSmps(smps_file);

# Parameter setting
DSPsolver.setLogLevel(1);                   # set print level
DSPsolver.setDdMasterNumCutsPerIter(10000); # set number of cuts per iteration (for the master problem)
DSPsolver.setWallLimit(wtime-60);           # set wallclock limit (sec)
DSPsolver.setDdMasterSolver(2);             # set solver type
DSPsolver.setDdAddFeasCuts(cut);             # enable feasibility cuts
DSPsolver.setDdAddOptCuts(cut);              # enable optimality cuts
if maxeval > 0
	DSPsolver.setDdMaxPrimsolEval(maxeval); # maximum number of solutions
else
	DSPsolver.setDdEvalUb(-1);
end
DSPsolver.setScipLimitsTime(300);

# solve problem using Dual Decomposition
DSPsolver.solve(DSP_SOLVER_DD)

# collect and print some results
if MPI.Comm_rank(MPI.COMM_WORLD) == 0

	# Write results
	f = open("output/$prefix.csv", "w");
	primal = DSPsolver.getPrimalBound();
	dual = DSPsolver.getDualBound();
	println(f, 
		"2,", 
		DSPsolver.getNumIterations(), ",", 
		primal, ",", 
		dual, ",", 
		DSPsolver.getSolutionTime(), ",",
		DSPsolver.getDdMasterTotalTime(), ",",
		DSPsolver.getDdLbTotalTime(), ",",
		DSPsolver.getDdUbTotalTime(), ",",
		DSPsolver.getDdCgTotalTime());
	close(f);

	# Print out simple results
	println("Solution status : ", DSPsolver.getSolutionStatus());
	println("Primal Bound    : ", DSPsolver.getPrimalBound())
	println("Dual Bound      : ", DSPsolver.getDualBound());
	println("CPU time (sec)  : ", DSPsolver.getDdCpuTime());

end

MPI.Finalize()
