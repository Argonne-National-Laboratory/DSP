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

using MPI, DSP

# initialize MPI
MPI.Init()

# read problem from SMPS
DSP.readSmps(smps_file);

# Parameter setting
DSP.setLogLevel(1);                   # set print level
DSP.setDdEvalUb(1);                   # evaluate upper bounds
DSP.setDdMasterNumCutsPerIter(10000); # set number of cuts per iteration (for the master problem)
DSP.setIterLimit(iterlim);            # set iteration limit
DSP.setWallLimit(wtime);              # set wallclock limit (sec)
DSP.setDdMasterSolver(solver_type);   # set solver type
DSP.setDdAddFeasCuts(fcut);           # enable feasibility cuts
DSP.setDdAddOptCuts(ocut);            # enable optimality cuts
DSP.setDdDualVarsLog(ddlog);          # enable logging dual variable values

# solve problem using Dual Decomposition
solve(DSP_SOLVER_DD)

# collect and print some results
if MPI.Comm_rank(MPI.COMM_WORLD) == 0

	# Write results
	f = open("output/$prefix.csv", "w");
	primal = DSP.getPrimalBound();
	dual = DSP.getDualBound();
	println(f, solver_type, ",", DSP.getNumIterations(), ",", primal, ",", dual, ",", DSP.getSolutionTime());
	close(f);

	if contains(smps_file, "sslp_15_45_15")
		writecsv("output/$prefix\_sub_objvals.csv", DSP.getDdSubproblemObjValues());
		writecsv("output/$prefix\_duals.csv", DSP.getDdDualBounds());
	end

	if ddlog == 1
		writecsv("output/$prefix\_dualvars.csv", DSP.getDdChangesOfMultiplier());
	end

	# Print out simple results
	println("Solution status : ", DSP.getSolutionStatus());
	println("Primal Bound    : ", DSP.getPrimalBound())
	println("Dual Bound      : ", DSP.getDualBound());
	println("CPU time (sec)  : ", DSP.getDdCpuTime());

end

MPI.Finalize()
