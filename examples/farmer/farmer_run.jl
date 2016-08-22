# Julia script for DSP with StochJuMP
# Kibaek Kim - ANL MCS 2014
# Farmer example from Birge and Louveaux book.

# model file
include("farmer_model.jl")

# Dsp solve types
solve_types = [:Dual, :Benders, :Extensive]

solve(m, solve_type = solve_types[1])

@show getprimobjval()
@show getdualobjval()
@show getprimvalue()
@show getdualvalue()
