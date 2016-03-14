# Julia script
# Kibaek Kim - ANL/MCS 2016

WORK_DIR  = "../examples/smps"
#INST_NAME = "sslp_15_45_5";
#INST_NAME = "sslp_10_50_50";
#INST_NAME = "farmer";
INST_NAME = "ieee118-2";
smps_file = "$WORK_DIR/$INST_NAME";

using MPI, DSPsolver

# initialize MPI
MPI.Init()

# retrieve MPI things
comm      = MPI.COMM_WORLD;
comm_rank = MPI.Comm_rank(comm);
comm_size = MPI.Comm_size(comm);

# read problem from SMPS
DSPsolver.readSmps(smps_file, true);

# set optionts
DSPsolver.setIntParam("LOG_LEVEL", 2);
DSPsolver.setIntParam("DD/MASTER_ALGO", REG_BUNDLE);
DSPsolver.setBoolParam("DD/ASYNC",true);
DSPsolver.setIntParam("ITER_LIM",5);

if comm_rank == 0
	tic();
end

# solve problem using extensive form
DSPsolver.solve(DSP_SOLVER_DD);

if comm_rank == 0
	elapsed_time = toq();
	@printf("solution time %f\n", elapsed_time);
end

for i=1:comm_size
	if i == 1 && comm_rank == 0
		# solution results
		dual = DSPsolver.getDualBound();
		# Print out simple results
		@printf("MPI rank %4d: Dual bound %+e\n", comm_rank, dual);
	elseif comm_rank == i - 1
		subprob_numiter  = DSPsolver.getDdNumSubproblemSolved();
		subprob_walltime = DSPsolver.getDdSubproblemWallTimes();
		subprob_totaltime = 0.0;
		@printf("MPI rank %4d: Number of iterations %d\n", comm_rank, subprob_numiter);
		for j=1:subprob_numiter
			@printf("  Iter %6d: solution time %.2f sec\n", j, subprob_walltime[j]);
			subprob_totaltime += subprob_walltime[j];
		end
		@printf("  Total solution time: %.2f sec\n", subprob_totaltime);
	end
	MPI.Barrier(comm);
end

MPI.Finalize();
