# Julia script for DSP with StochJuMP
# Kibaek Kim - ANL MCS 2014
# Farmer example from Birge and Louveaux book.

using MPI, StochJuMP, DSPsolver

if ARGS[1] == "DD"
	MPI.Init();
end

# data file
include("farmer_data.jl")

# model file
include("farmer_model.jl")

# load problem to model object
DSPsolver.loadProblem(m);

# set parameters
DSPsolver.setIntParam("LOG_LEVEL",1);
DSPsolver.setIntParam("DD/ITER_LIM",100);
DSPsolver.setDblParam("SCIP/GAP_TOL",0.0);

# solve problem
if ARGS[1] == "DE"
	DSPsolver.solve(DSP_SOLVER_DE);
elseif ARGS[1] == "BD"
	DSPsolver.setIntParam("BD/NUM_CORES",1);
	DSPsolver.solve(DSP_SOLVER_BD);
elseif ARGS[1] == "DD"
	DSPsolver.setIntParam("DD/FEAS_CUTS",1);
	DSPsolver.setIntParam("DD/OPT_CUTS",1);
	DSPsolver.setIntParam("DD/EVAL_UB",1);
	DSPsolver.setIntParam("DD/MASTER_ALGO",IPM_FEAS);
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

