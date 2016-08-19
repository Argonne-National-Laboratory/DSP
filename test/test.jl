# Julia script for testing DSP
# Kibaek Kim - 2016

using Dsp

# parameter file
myparam = joinpath(dirname(@__FILE__),"params.set")

# smps directory
smps_dir = joinpath(dirname(@__FILE__), "../examples/smps")

# model file
include("../examples/farmer/farmer_model.jl")

# Dsp solve types
solve_types = [:Dual, :Benders, :Extensive]

for st in solve_types
    println("\n################################################################################")
    println("# Solve farmer using ", st)
    println("################################################################################\n")
    solve(m, solve_type = st, param = myparam)
    @show getprimobjval()
    @show getdualobjval()
    @show Dsp.model.colVal
    @show Dsp.model.rowVal
    @show getsolutiontime()
end

smps_files = readdir(smps_dir)
for f in smps_files
    if contains(f, ".cor") == true
        pos = search(f, '.') - 1
        fname = joinpath(smps_dir, f[1:pos])
        readSmps(fname)
        for st in solve_types
            println("\n################################################################################")
            println("# Solve $fname using ", st)
            println("################################################################################\n")
            optimize(solve_type = st, param = myparam)
            @show getprimobjval()
            @show getdualobjval()
            @show getsolutiontime()
        end
    end
end