# Julia script
# Kibaek Kim - ANL MCS 2016

using MPI, Dsp

if length(ARGS) < 1
	error("Please give an SMPS file name.")
end

smps = ARGS[1]

# initialize MPI
MPI.Init()

readSmps(smps)
solve()

MPI.Finalize()
