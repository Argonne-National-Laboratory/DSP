# Julia script for testing DSP
# Kibaek Kim - 2016

using MPI, Dsp

MPI.Init()

# parameter file
myparam = joinpath(dirname(@__FILE__),"params.set")

# smps directory
smps_dir = joinpath(dirname(@__FILE__), "../examples/smps")

# model file
include("../examples/farmer/farmer_model.jl")

# Dsp solve types
solve_types = [:Dual, :Benders]

# for st in solve_types
#     println("\n################################################################################")
#     println("# Solve farmer using ", st)
#     println("################################################################################\n")
#     solve(m, solve_type = st, param = myparam, comm = MPI.COMM_WORLD)
#     if MPI.Comm_rank(MPI.COMM_WORLD) == 0
#         @show getprimobjval()
#         @show getdualobjval()
#         @show getsolutiontime()
#     end
# end

smps_files = readdir(smps_dir)
# for f in smps_files
f = smps_files[1]
    if contains(f, ".cor") == true
        pos = search(f, '.') - 1
        fname = joinpath(smps_dir, f[1:pos])
        readSmps(fname)
        for st in solve_types
            println("\n################################################################################")
            println("# Solve $fname using ", st)
            println("################################################################################\n")
            optimize(solve_type = st, param = myparam, comm = MPI.COMM_WORLD)
            if MPI.Comm_rank(MPI.COMM_WORLD) == 0
                @show getprimobjval()
                @show getdualobjval()
                @show getsolutiontime()
            end
        end
    end
# end

MPI.Finalize()