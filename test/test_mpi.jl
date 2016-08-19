# Julia script for testing DSP
# Kibaek Kim - 2016

using MPI, Dsp

test_name    = AbstractString[]
test_method  = Symbol[]
test_status  = Symbol[]
test_primobj = Float64[]
test_dualobj = Float64[]
test_time    = Float64[]

MPI.Init()

# parameter file
myparam = joinpath(dirname(@__FILE__),"params.set")

# smps directory
smps_dir = joinpath(dirname(@__FILE__), "../examples/smps")

# model file
include("../examples/farmer/farmer_model.jl")

# Dsp solve types
solve_types = [:Dual, :Benders]

for st in solve_types
    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        println("\n################################################################################")
        println("# Solve farmer using ", st)
        println("################################################################################\n")
    end
    status = solve(m, solve_type = st, param = myparam, comm = MPI.COMM_WORLD)
    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        @show getprimobjval()
        @show getdualobjval()
        @show getsolutiontime()
        push!(test_name, "farmer(JuMP)")
        push!(test_method, st)
        push!(test_status, status)
        push!(test_primobj, getprimobjval())
        push!(test_dualobj, getdualobjval())
        push!(test_time, getsolutiontime())
    end
end

smps_files = readdir(smps_dir)
for f in smps_files
    if contains(f, ".cor") == true
        pos = search(f, '.') - 1
        fname = joinpath(smps_dir, f[1:pos])
        readSmps(fname)
        for st in solve_types
            if MPI.Comm_rank(MPI.COMM_WORLD) == 0
                println("\n################################################################################")
                println("# Solve $fname using ", st)
                println("################################################################################\n")
            end
            status = optimize(solve_type = st, param = myparam, comm = MPI.COMM_WORLD)
            if MPI.Comm_rank(MPI.COMM_WORLD) == 0
                @show getprimobjval()
                @show getdualobjval()
                @show getsolutiontime()
                push!(test_name, f[1:pos])
                push!(test_method, st)
                push!(test_status, status)
                push!(test_primobj, getprimobjval())
                push!(test_dualobj, getdualobjval())
                push!(test_time, getsolutiontime())
            end
        end
    end
end

if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    for i in 1:length(test_name)
        @printf("%15s  %10s  %15s  %+10e  %+10e  %6.2f\n", 
            test_name[i], test_method[i], test_status[i], test_primobj[i], test_dualobj[i], test_time[i])
    end
end

MPI.Finalize()