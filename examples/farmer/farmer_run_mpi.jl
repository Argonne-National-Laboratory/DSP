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
DSPsolver.setIntParam("DD/ITER_LIM",100);
DSPsolver.setDblParam("SCIP/GAP_TOL",0.0);

init_solution = [170.0,79.99999999999997,250.00000000000003,0.0,0.0,310.0,47.999999999999886,6000.0,9.094947017729282e-13,0.0,0.0,225.0,-8.526512829121202e-14,5000.000000000001,0.0,0.0,48.000000000000085,140.0,0.0,4000.0000000000005,0.0];
DSPsolver.setSolution(init_solution);

# solve problem
DSPsolver.solve(DSP_SOLVER_BD_MPI);
#DSPsolver.solveBdMpi();

if MPI.Comm_rank(MPI.COMM_WORLD) == 0
	# print some results
	println("Solution Status: ", DSPsolver.getSolutionStatus());
	println("Primal Bound   : ", DSPsolver.getPrimalBound());
	println("Dual Bound     : ", DSPsolver.getDualBound());

	print(DSPsolver.getSolution());
end

MPI.Finalize();

