# Julia script for DSP with StochJuMP
# Kibaek Kim - ANL MCS 2014
# Farmer example from Birge and Louveaux book.

# model file
include("farmer_model.jl")

# Dsp solve types
solve_types = [:Dual, :Benders, :Extensive]

solve(m, solve_type = solve_types[2])

@show getprimobjval() # Dsp.model.primVal
@show getdualobjval() # Dsp.model.dualVal
@show Dsp.model.colVal
@show Dsp.model.rowVal