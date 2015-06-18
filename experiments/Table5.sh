#!/bin/bash

# algorithmic parameters
nscenarios=(4 8 16 32 64)

for i in {0..4}
do
	# executing script
	exec_prefix="suc${nscenarios[i]}-de"
	exec_script="julia suc_run_de.jl ${nscenarios[i]} ${exec_prefix}"
	mpirun -n 1 ${exec_script}
done
