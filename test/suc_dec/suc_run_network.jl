println("-----------------------------------------------------");
println("Julia script for Unit Commitment Problem.");
println("Kibaek Kim - 2015 ANL MCS");
println("Note:");
println(" The original problem was written in AMPL scripts,");
println(" which were provided by Changhyeok Lee and Cong Liu.");
println("-----------------------------------------------------");

if length(ARGS) < 1
	println("Arguments must be specified");
	exit();
end

using MPI, StochJuMP, DSPsolver

# Algorithmic parameters
nScenarios = 1;                  # number of scenarios
# nScenarios = 4;                  # number of scenarios
nScenariosToCheck = 1;           # number of scenarios to check for feasibility for each scenario subproblem
scenariosToSolve = 1:nScenarios; # set of scenarios solved by each processor
scenariosToRead = 1:nScenarios;  # set of scenarios read by each processor

# MPI
comm_size = -1;
comm_rank = -1;

if ARGS[1] == "DD"
	# Initialize MPI
	MPI.Init();

	# Process information
	comm_size = MPI.Comm_size(MPI.COMM_WORLD);
	comm_rank = MPI.Comm_rank(MPI.COMM_WORLD);

	# create scenario indices to solve and read
	scenariosToSolve = Int[];
	scenariosToRead = Int[];
	for s in (comm_rank+1):comm_size:nScenarios
		push!(scenariosToSolve, s);
		push!(scenariosToRead, s);
	end
#	for s in scenariosToSolve
#		for i in 1:nScenariosToCheck
#			if s + i > nScenarios
#				push!(scenariosToRead, s + i - nScenarios)
#			else
#				push!(scenariosToRead, s + i)
#			end
#		end
#	end
	scenariosToRead = sort(unique(scenariosToRead));
	for r in 1:comm_size
		if r - 1 == comm_rank
			println("Processor ", comm_rank);
			println("  to solve ", scenariosToSolve);
			println("  to read  ", scenariosToRead);
		end
		MPI.Barrier(MPI.COMM_WORLD);
	end
end

# --------------
# Initialize DSP
# --------------

if length(ARGS) >= 2
	nclusters = int(ARGS[2]);
else
	nclusters = 3;
end

# Unit commitment model
print("Reading model file, ")
@time include("suc_mod_dec_network.jl")

# garbage collecting
gc()

DSPsolver.setLogLevel(4)
DSPsolver.setNumCores(1)

if ARGS[1] == "DE"
	DSPsolver.solveDe()
	println("Solution status: ", DSPsolver.getSolutionStatus())
	println("Primal Bound   : ", DSPsolver.getPrimalBound())
	println("Dual Bound     : ", DSPsolver.getDualBound())

	println(DSPsolver.getSolution())
elseif ARGS[1] == "DD"
	# DSPsolver.setDdMasterSolver(3) # Subgradient
	DSPsolver.setDdMasterSolver(2)
	DSPsolver.setDdMasterNumCutsPerIter(0)
	DSPsolver.setDdAddFeasCuts(-1)
	DSPsolver.setDdAddOptCuts(-1)
	DSPsolver.setDdEvalUb(-1)
	DSPsolver.setDdDualVarsLog( 1)
	DSPsolver.setDdTrustRegionSize(1000)
	DSPsolver.setDdDisableTrustRegionDecrease(true)
	DSPsolver.solveDd(MPI.COMM_WORLD)

	if MPI.Comm_rank(MPI.COMM_WORLD) == 0
		writecsv("dd_itertime.csv", DSPsolver.getDdIterTime())
		writecsv("dd_subtime.csv", DSPsolver.getDdSubprobTime())
		writecsv("dd_master_objvals.csv", DSPsolver.getDdMasterObjValues())
		writecsv("dd_sub_objvals.csv", DSPsolver.getDdSubproblemObjValues())
		writecsv("dd_primals.csv", DSPsolver.getDdPrimalBounds())
		writecsv("dd_duals.csv", DSPsolver.getDdDualBounds())
		println("Solution status: ", DSPsolver.getSolutionStatus())
		println("Primal Bound   : ", DSPsolver.getPrimalBound())
		println("Dual Bound     : ", DSPsolver.getDualBound())
	end
end

if ARGS[1] == "DD"
	# Initialize MPI
	MPI.Finalize();
end
