# Julia script
# Kibaek Kim - ANL MCS 2016

using Dsp

if length(ARGS) < 1
	error("Please give an SMPS file name.")
end

smps = ARGS[1]

readSmps(smps)

# type of solution methods
solve_types = [:Dual, :Benders, :Extensive]

# Default parameter file
myparam = joinpath(dirname(@__FILE__),"../parameters/default.txt")

optimize(solve_type = solve_types[1], param = myparam)

@show getprimobjval() # Dsp.model.primVal
@show getdualobjval() # Dsp.model.dualVal
@show Dsp.model.colVal
@show Dsp.model.rowVal