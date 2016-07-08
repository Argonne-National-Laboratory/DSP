import PyDsp
from mpi4py import MPI

# create DSP environment
dsp = PyDsp.createEnv()

# read model from SMPS files
PyDsp.readSmps(dsp, "../examples/smps/farmer")

# MPI information
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

# solve Benders decomposition
PyDsp.solveBd(dsp, 1, comm)
comm.Barrier()

# free DSP environment
PyDsp.freeEnv(dsp)
