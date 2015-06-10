# Julia interface for DSP + MPI
# Kibaek Kim - 2015 ANL MCS

include("Dsp.jl")

# Check if MPI is installed.
if Pkg.installed("MPI") == nothing
	Pkg.add("MPI")
end
using MPI

function solveDd(env::DSP, comm)
	@dsp_ccall("solveDd", Void, (Ptr{Void}, Cint), env.p, convert(Cint, comm.val))
end

