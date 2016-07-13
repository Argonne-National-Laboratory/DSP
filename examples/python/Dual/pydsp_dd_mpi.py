import argparse
import PyDsp
from mpi4py import MPI

__author__ = 'Kibaek Kim'
parser = argparse.ArgumentParser(description='This is a PyDsp example for dual decomposition.')
parser.add_argument('-s', '--smps', help='SMPS file name', required=True)
parser.add_argument('-p', '--param', help='Parameter file name', required=False)
args = parser.parse_args()

# create DSP environment
dsp = PyDsp.createEnv()

# read model from SMPS files
PyDsp.readSmps(dsp, args.smps)

# read parameter file
if args.param is not None:
	PyDsp.readParamFile(dsp, args.param)

# MPI information
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

# solve dual decomposition
PyDsp.solveDdMpi(dsp, comm)
comm.Barrier()

# solution status
if rank == 0:
	status = PyDsp.getStatus(dsp)
	print ("Solution status: %d" % status)
	
	# collect solutions
	ncols = PyDsp.getTotalNumCols(dsp)
	print ("Number of variables: %d" % ncols)
	print ("Primal bound: %+e" % PyDsp.getPrimalBound(dsp))
	print ("Dual bound  : %+e" % PyDsp.getDualBound(dsp))
	'''
	x_val = PyDsp.doubleArray(ncols)
	PyDsp.getSolution(dsp, ncols, x_val)
	for i in xrange(ncols):
		if abs(x_val[i]) > 1.0e-6:
			print ("x[%d]\t%+e" % (i,x_val[i]))
	'''

# free DSP environment
PyDsp.freeEnv(dsp)
