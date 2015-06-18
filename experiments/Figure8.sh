#!/bin/bash

# algorithmic parameters
nscenarios=8
solver=(3 2)
cuts=(-1 -1)

# mpirun setting
nmpinodes=8

for i in {0..1}
do
	# executing script
	exec_prefix="suc${nscenarios}-${solver[i]}${cuts[i]}"
	exec_script="julia suc_run.jl ${nscenarios} ${solver[i]} ${cuts[i]} ${exec_prefix}"

	mpirun -n ${nmpinodes} ${exec_script}
done
