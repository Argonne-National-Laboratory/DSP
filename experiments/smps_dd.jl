# Julia script of reproducing results for Table 2, Figures 1 and 2 in MPC paper.
# Kibaek Kim - ANL MCS 2015

# Arguments from command line
smps_file   = ARGS[1];      # SMPS file name
solver_type = int(ARGS[2]); # 0 = Simplex, 2 = Primal-dual interior point, 3 = Subgradient
fcut        = int(ARGS[3]); # add feasibility cut: -1 = No
ocut        = int(ARGS[4]); # add optimality cut : -1 = No
ddlog       = int(ARGS[5]); # Log dual variable values?
iterlim     = int(ARGS[6]); # iteration limit
wtime       = int(ARGS[7]); # wall time limit
prefix      = ARGS[8];      # output prefix

using MPI, DSPsolver

# initialize MPI
MPI.Init()

# read problem from SMPS
DSPsolver.readSmps(smps_file);

# Parameter setting
DSPsolver.setLogLevel(1);                   # set print level
DSPsolver.setDdMasterNumCutsPerIter(10000); # set number of cuts per iteration (for the master problem)
DSPsolver.setIterLimit(iterlim);            # set iteration limit
DSPsolver.setWallLimit(wtime);              # set wallclock limit (sec)
DSPsolver.setDdMasterSolver(solver_type);   # set solver type
DSPsolver.setDdAddFeasCuts(fcut);           # enable feasibility cuts
DSPsolver.setDdAddOptCuts(ocut);            # enable optimality cuts
DSPsolver.setDdDualVarsLog(ddlog);          # enable logging dual variable values

# solve problem using Dual Decomposition
DSPsolver.solve(DSP_SOLVER_DD)

# collect and print some results
if MPI.Comm_rank(MPI.COMM_WORLD) == 0

	# Write results
	f = open("output/$prefix.csv", "w");
	primal = DSPsolver.getPrimalBound();
	dual = DSPsolver.getDualBound();
	println(f, solver_type, ",", DSPsolver.getNumIterations(), ",", primal, ",", dual, ",", DSPsolver.getSolutionTime());
	close(f);

	if contains(smps_file, "sslp_15_45_15")
		writecsv("output/$prefix\_sub_objvals.csv", DSPsolver.getDdSubproblemObjValues());
		writecsv("output/$prefix\_duals.csv", DSPsolver.getDdDualBounds());
	end

	if ddlog == 1
		writecsv("output/$prefix\_dualvars.csv", DSPsolver.getDdChangesOfMultiplier());
	end

	# Print out simple results
	println("Solution status : ", DSPsolver.getSolutionStatus());
	println("Primal Bound    : ", DSPsolver.getPrimalBound())
	println("Dual Bound      : ", DSPsolver.getDualBound());
	println("CPU time (sec)  : ", DSPsolver.getDdCpuTime());

end

MPI.Finalize()
