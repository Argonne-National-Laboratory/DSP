#!/bin/bash

# algorithmic parameters
nscenarios=(4 8 16 32 64)

# mpirun setting
nmpinodes=(4 8 16 32 64)
npernode=(4 8 16 16 16)

# cluster setting
nnodes=(1 1 1 2 4)
ncores=(4 8 16 16 16)

for i in {0..4}
do
	# executing script
	exec_prefix="suc${nscenarios[i]}_fscip"
	exec_args="mps/suc${nscenarios[i]}.mps -fsol output/suc${nscenarios[i]}.fscip.sol"
	exec_script="./fscip default.set ${exec_args}"

	cp pbs/template.pbs pbs/${exec_prefix}.pbs
	echo "mpirun -machinefile hostfile -n ${nmpinodes[i]} -ppn ${npernode[i]} ${exec_script}" >> pbs/${exec_prefix}.pbs
	qsub -N ${exec_prefix} -l nodes=${nnodes[i]}:ppn=${ncores[i]},walltime=6:30:00 pbs/${exec_prefix}.pbs
done
