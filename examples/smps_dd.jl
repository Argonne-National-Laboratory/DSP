# Julia script for DSP (Dual Decomposition)
# Kibaek Kim - ANL MCS 2014

# DSP solver
include("../julia/Dsp-MPI.jl")

# create DSP environment
dsp = DSP()

# read problem from SMPS
#readSmps(dsp, "./smps/dcap233_200")
readSmps(dsp, "./smps/sslp_5_25_50")

setLogLevel(dsp, 2);                 # set print level
#setIterLimit(dsp, 100);             # set iteration limit
#setWallLimit(dsp, 1000000);          # set wallclock limit (sec)
setDdMasterNumCutsPerIter(dsp, 100); # number of cuts per iteration
setDdMasterSolver(dsp, 2);           # use IPM
setDdAddFeasCuts(dsp,-1);             # add feasibility cuts
setDdAddOptCuts(dsp,-1);              # add optimality cuts
setDdEvalUb(dsp,-1);                  # evaluate upper bounds

# initialize MPI
MPI.Init()

# solve problem using Dual Decomposition
solveDd(dsp, MPI.COMM_WORLD)

# collect and print some results
if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    println("Solution status : ", getSolutionStatus(dsp))
    println("Primal Bound    : ", getPrimalBound(dsp))
    println("Dual Bound      : ", getDualBound(dsp))
    println("CPU time (sec)  : ", getDdCpuTime(dsp))
end

MPI.Finalize()
