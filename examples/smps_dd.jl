# Julia script of reproducing results for Table 2, Figures 1 and 2 in MPC paper.
# Kibaek Kim - ANL MCS 2015

using MPI, DSPsolver

smps_file   = "/home/kibaekkim/DSP/examples/smps/sslp_5_25_50";
solver_type = IPM_FEAS;
fcut        = 1;
ocut        = 1;
ddlog       = false;
iterlim     = 100;
wtime       = 100000;
prefix      = "sslp_5_25_50";

# Arguments from command line
if length(ARGS) == 8
	smps_file   = ARGS[1];      # SMPS file name
	solver_type = parse(Int,ARGS[2]); # 0 = Simplex, 1 = IPM, 2=IPM-FEAS, 3 = DSBM, 4 = Subgradient
	fcut        = parse(Int,ARGS[3]); # add feasibility cut: 0 = No
	ocut        = parse(Int,ARGS[4]); # add optimality cut : 0 = No
	ddlog       = parse(Int,ARGS[5]); # Log dual variable values?
	iterlim     = parse(Cuchar,ARGS[6]); # iteration limit
	wtime       = parse(Int,ARGS[7]); # wall time limit
	prefix      = ARGS[8];      # output prefix
end

# initialize MPI
MPI.Init()

# read problem from SMPS
DSPsolver.readSmps(smps_file);

# Parameter setting
DSPsolver.setIntParam("LOG_LEVEL",2);
DSPsolver.setIntParam("DD/NUM_CUTS_PER_ITER",10000);
DSPsolver.setIntParam("ITER_LIM",iterlim);
DSPsolver.setDblParam("WALL_LIM",wtime);
DSPsolver.setIntParam("DD/MASTER_ALGO",solver_type);
DSPsolver.setIntParam("DD/FEAS_CUTS",fcut);
DSPsolver.setIntParam("DD/OPT_CUTS",ocut);
DSPsolver.setBoolParam("DD/LOG_DUAL_VARS",ddlog);
DSPsolver.setDblParam("SCIP/GAP_TOL",100.0);

# solve problem using Dual Decomposition
DSPsolver.solve(DSP_SOLVER_DD)

# collect and print some results
if MPI.Comm_rank(MPI.COMM_WORLD) == 0

	# Write results
	f = open("$prefix.csv", "w");
	primal = DSPsolver.getPrimalBound();
	dual = DSPsolver.getDualBound();
	println(f, solver_type, ",", DSPsolver.getNumIterations(), ",", primal, ",", dual, ",", DSPsolver.getSolutionTime());
	close(f);

	if contains(smps_file, "sslp_15_45_15")
		writecsv("$prefix\_sub_objvals.csv", DSPsolver.getDdSubproblemObjValues());
		writecsv("$prefix\_duals.csv", DSPsolver.getDdDualBounds());
	end

	if ddlog == 1
		writecsv("$prefix\_dualvars.csv", DSPsolver.getDdChangesOfMultiplier());
	end

	# Print out simple results
	println("Solution status : ", DSPsolver.getSolutionStatus());
	println("Primal Bound    : ", DSPsolver.getPrimalBound())
	println("Dual Bound      : ", DSPsolver.getDualBound());
	println("CPU time (sec)  : ", DSPsolver.getDdCpuTime());

end

MPI.Finalize()
