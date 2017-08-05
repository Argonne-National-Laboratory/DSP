# Julia script for DSP with StochJuMP
# Kibaek Kim - ANL MCS 2014
# Farmer example from Birge and Louveaux book.

# model file
include("farmer_model.jl")

# Dsp solve types
solve_types = [:Dual, :Benders, :Extensive]

# Default parameter file
myparam = joinpath(dirname(@__FILE__),"../../parameters/default.txt")

solve(m, solve_type = solve_types[1], param = myparam)

@show getprimobjval()
@show getdualobjval()
@show getprimvalue()
@show getdualvalue()
