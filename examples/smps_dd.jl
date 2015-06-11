# Julia script for DSP (Dual Decomposition)
# Kibaek Kim - ANL MCS 2014

using DSP, MPI

MPI.Init() # initialize MPI

DSP.readSmps("./smps/dcap233_200"); # read SMPS files
DSP.setLogLevel(1);                 # set print level
DSP.solve(DSP_SOLVER_DD);           # solve with DD

# collect and print some results
if MPI.Comm_rank(MPI.COMM_WORLD) == 0
	println("Solution Time   : ", DSP.getSolutionTime());
	println("Solution status : ", DSP.getSolutionStatus());
	println("Primal Bound    : ", DSP.getPrimalBound());
	println("Dual Bound      : ", DSP.getDualBound());
	println("Iterations      : ", DSP.getNumIterations());
end

MPI.Finalize(); # finalize MPI
