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

status = solve(m, solve_type = solve_types[2])

if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    @show getprimobjval()
    @show getdualobjval()
    @show getprimvalue()
    @show getdualvalue()
end

for i = 1:MPI.Comm_size(MPI.COMM_WORLD)
    if MPI.Comm_rank(MPI.COMM_WORLD) == i - 1
        @show MPI.Comm_rank(MPI.COMM_WORLD)
        for b in m.ext[:DspBlocks].children
            @show b.first
            @show b.second.colVal
        end
    end
    MPI.Barrier(MPI.COMM_WORLD)
end
# finalize
MPI.Finalize()
