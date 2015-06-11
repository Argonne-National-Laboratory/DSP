# Julia script for DSP with StochJuMP
# Kibaek Kim - ANL MCS 2014
# Farmer example from Birge and Louveaux book.

using MPI, StochJuMP, DSP

if ARGS[1] == "DD"
	MPI.Init();
end

include("farmer_data.jl");  # data file
include("farmer_model.jl"); # model file

DSP.loadProblem(m); # load problem to model object
DSP.setLogLevel(1); # set print level


# solve problem
if ARGS[1] == "DE"
	DSP.solve(DSP_SOLVER_DE);
elseif ARGS[1] == "BD"
	DSP.solve(DSP_SOLVER_BD);
elseif ARGS[1] == "DD"
	DSP.solve(DSP_SOLVER_DD);
end

if MPI.Initialized() == false || MPI.Comm_rank(MPI.COMM_WORLD) == 0
	println("Solution Time   : ", DSP.getSolutionTime());
	println("Solution status : ", DSP.getSolutionStatus());
	println("Primal Bound    : ", DSP.getPrimalBound());
	println("Dual Bound      : ", DSP.getDualBound());
	println("Iterations      : ", DSP.getNumIterations());
end

if ARGS[1] == "DD"
	MPI.Finalize();
end

