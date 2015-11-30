# Julia script for DSP with StochJuMP
# Kibaek Kim - ANL MCS 2014
# Farmer example from Birge and Louveaux book.

using MPI, StochJuMP, DSPsolver

MPI.Init();

# data file
include("farmer_data.jl")

# model file
include("farmer_model.jl")

# load problem to model object
DSPsolver.loadProblem(m);

# set parameters
DSPsolver.setIntParam("LOG_LEVEL",3);
DSPsolver.setIntParam("ITER_LIM",100);
DSPsolver.setDblParam("SCIP/GAP_TOL",0.0);

# solve problem
# DSPsolver.solve(DSP_SOLVER_BD_MPI);
DSPsolver.solveBdMpi(3);

if MPI.Comm_rank(MPI.COMM_WORLD) == 0
	# print some results
	println("Solution Status: ", DSPsolver.getSolutionStatus());
	println("Primal Bound   : ", DSPsolver.getPrimalBound());
	println("Dual Bound     : ", DSPsolver.getDualBound());
end

MPI.Finalize();

