# Julia script
# Kibaek Kim - ANL MCS 2015
# updated on 2016

using MPI, Dsp

# Initialize MPI
MPI.Init()

# Algorithmic parameters
nScenarios  = 2 # number of scenarios
if length(ARGS) == 1
	nScenarios  = int(ARGS[1])
end

# Unit commitment model
include("suc_mod.jl")

solve(m, solve_type = :Dual, comm = MPI.COMM_WORLD)

if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    @show getprimobjval() # Dsp.model.primVal
    @show getdualobjval() # Dsp.model.dualVal
    @show Dsp.model.rowVal
end

MPI.Finalize()

