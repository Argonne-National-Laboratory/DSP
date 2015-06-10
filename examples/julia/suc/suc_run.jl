println("-----------------------------------------------------");
println("Julia script for Unit Commitment Problem.");
println("Kibaek Kim - 2015 ANL MCS");
println("Note:");
println(" The original problem was written in AMPL scripts,");
println(" which were provided by Changhyeok Lee and Cong Liu.");
println("-----------------------------------------------------");

# StochJuMP
include("../../../julia/StochJuMP.jl")

# DSP
include("../../../julia/Dsp-MPI.jl")

# Algorithmic parameters
nScenarios = 4;                  # number of scenarios
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
			println("Procesor ", comm_rank);
			println("  to solve ", scenariosToSolve);
			println("  to read  ", scenariosToRead);
		end
		MPI.Barrier(MPI.COMM_WORLD);
	end
end

# --------------
# Initialize DSP
# --------------
dsp = DSP()

# Unit commitment model
print("Reading model file, ")
@time include("suc_mod_with_fastgen.jl")
#@time include("suc_mod.jl")

# garbage collecting
gc()

setLogLevel(dsp, 4)
setNumCores(dsp, 1)

if ARGS[1] == "DE"
	solveDe(dsp)
	println("Solution status: ", getSolutionStatus(dsp))
	println("Primal Bound   : ", getPrimalBound(dsp))
	println("Dual Bound     : ", getDualBound(dsp))
elseif ARGS[1] == "BD"
	# augmented scenarios
	setIntRelax(dsp,1);
	augs = [0,1,2];
	naugs = length(augs);
	setBdAugScenarios(dsp, naugs, augs)
	solveBd(dsp, 1)
	println("Solution status: ", getSolutionStatus(dsp))
	println("Primal Bound   : ", getPrimalBound(dsp))
	println("Dual Bound     : ", getDualBound(dsp))
elseif ARGS[1] == "DD"
	setDdMasterSolver(dsp, 0)
	setDdMasterNumCutsPerIter(dsp,1)
	setDdAddFeasCuts(dsp,1)
	setDdAddOptCuts(dsp,1)
	setDdEvalUb(dsp,1)
	#setIterLimit(dsp,10)
	solveDd(dsp, MPI.COMM_WORLD)
	if MPI.Comm_rank(MPI.COMM_WORLD) == 0
		writecsv("dd_itertime.csv", getDdIterTime(dsp))
	        writecsv("dd_subtime.csv", getDdSubprobTime(dsp))
	        writecsv("dd_master_objvals.csv", getDdMasterObjValues(dsp))
        	writecsv("dd_sub_objvals.csv", getDdSubproblemObjValues(dsp))
	        writecsv("dd_primals.csv", getDdPrimalBounds(dsp))
        	writecsv("dd_duals.csv", getDdDualBounds(dsp))
		println("Solution status: ", getSolutionStatus(dsp))
		println("Primal Bound   : ", getPrimalBound(dsp))
		println("Dual Bound     : ", getDualBound(dsp))
	end
end

if ARGS[1] == "DD"
	# Initialize MPI
	MPI.Finalize();
end
