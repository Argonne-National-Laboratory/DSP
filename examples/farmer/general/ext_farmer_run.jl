# Julia script for DSP with StochJuMP
# Kibaek Kim - ANL MCS 2014
# Farmer example from Birge and Louveaux book.

using MPI, StochJuMP, DSPsolver

if ARGS[1] == "DD"
	MPI.Init();
end

# data file
include("../farmer_data.jl")

# model file
include("ext_farmer_model.jl")

# load problem to model object
DSPsolver.loadProblem(m);

# load decomposition
DSPsolver.loadDecomposition(m, cc);

# set print level
DSPsolver.setLogLevel(1);

# solve problem
if ARGS[1] == "DE"
	DSPsolver.solve(DSP_SOLVER_DE);
elseif ARGS[1] == "BD"
	DSPsolver.setNumCores(1); # set number of cores
	DSPsolver.solve(DSP_SOLVER_BD);
elseif ARGS[1] == "DD"
	DSPsolver.setDdAddFeasCuts(-1);
	DSPsolver.setDdAddOptCuts(-1);
	DSPsolver.setDdEvalUb(-1);
	DSPsolver.solve(DSP_SOLVER_DD);
end

if MPI.Initialized() == false || MPI.Comm_rank(MPI.COMM_WORLD) == 0
	# print some results
	println("Solution Status: ", DSPsolver.getSolutionStatus());
	println("Primal Bound   : ", DSPsolver.getPrimalBound());
	println("Dual Bound     : ", DSPsolver.getDualBound());
end

if ARGS[1] == "DD"
	MPI.Finalize();
end

