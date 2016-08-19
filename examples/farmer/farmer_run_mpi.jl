# Julia script for DSP with StochJuMP
# Kibaek Kim - ANL MCS 2014
# Farmer example from Birge and Louveaux book.

using MPI

# initialize
MPI.Init();

# model file
include("farmer_model.jl")

# Dsp solve types
solve_types = [:Dual, :Benders]

status = solve(m, solve_type = solve_types[2], comm = MPI.COMM_WORLD)

if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    @show getprimobjval() # Dsp.model.primVal
    @show getdualobjval() # Dsp.model.dualVal
    @show Dsp.model.colVal
    @show Dsp.model.rowVal
end

# finalize
MPI.Finalize()
