#!/bin/bash

# algorithmic parameters
nscenarios=8
solver=(3 2)
cuts=(-1 -1)

# mpirun setting
nmpinodes=8
npernode=8

# cluster setting
nnodes=1
ncores=8

for i in {0..1}
do
	# executing script
	exec_prefix="suc${nscenarios}-${solver[i]}${cuts[i]}"
	exec_script="julia suc_run.jl ${nscenarios} ${solver[i]} ${cuts[i]} ${exec_prefix}"

	cp pbs/template.pbs pbs/${exec_prefix}.pbs
	echo "mpirun -machinefile hostfile -n ${nmpinodes} -ppn ${npernode} ${exec_script}" >> pbs/${exec_prefix}.pbs
	qsub -N ${exec_prefix} -l nodes=${nnodes}:ppn=${ncores},walltime=6:30:00 pbs/${exec_prefix}.pbs
done
