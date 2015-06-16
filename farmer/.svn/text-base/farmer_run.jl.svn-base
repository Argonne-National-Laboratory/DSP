# Julia script for DSP with StochJuMP
# Kibaek Kim - ANL MCS 2014
# Farmer example from Birge and Louveaux book.

# read include path from environment variable
DSP_INC = ENV["DSP_INC"];

# modeling language, StochJuMP
include("$DSP_INC/StochJuMP.jl")

# solver, DSP
include("$DSP_INC/Dsp-MPI.jl")

if ARGS[1] == "DD"
	MPI.Init();
end

# data file
include("farmer_data.jl")

# model file
include("farmer_model.jl")

# create DSP environment
dsp = DSP();

# load problem to model object
loadProblem(dsp, m);

# set print level
setLogLevel(dsp, 1);

# solve problem
if ARGS[1] == "DE"
	solveDe(dsp);
elseif ARGS[1] == "BD"
	setNumCores(dsp, 3);
	solveBd(dsp, 1);
elseif ARGS[1] == "DD"
	setDdMasterSolver(dsp, 2);
	setDdMasterNumCutsPerIter(dsp, 1);
	setDdAddFeasCuts(dsp,1);
	setDdAddOptCuts(dsp,1);
	setDdEvalUb(dsp,1);
	solveDd(dsp, MPI.COMM_WORLD);
end

if MPI.Initialized() == false || MPI.Comm_rank(MPI.COMM_WORLD) == 0
	# print some results
	println("Solution Status: ", getSolutionStatus(dsp))
	println("Primal Bound   : ", getPrimalBound(dsp))
	println("Dual Bound     : ", getDualBound(dsp))
end

if ARGS[1] == "DD"
	MPI.Finalize();
end

