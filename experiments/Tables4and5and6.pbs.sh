#!/bin/bash

# algorithmic parameters
nscenarios=(4 8 16 32 64)
solver=(0 0 2 2 3)
cuts=(-1 1 -1 1 -1)

# mpirun setting
nmpinodes=(4 8 16 32 64)
npernode=(4 8 16 16 16)

# cluster setting
nnodes=(1 1 1 2 4)
ncores=(4 8 16 16 16)

for i in {0..4}
do
	for j in {0..4}
	do
		# executing script
		exec_prefix="suc${nscenarios[i]}_${solver[j]}${cuts[j]}"
		exec_script="julia suc_run.jl ${nscenarios[i]} ${solver[j]} ${cuts[j]} ${exec_prefix}"

		cp pbs/template.pbs pbs/${exec_prefix}.pbs
		echo "mpirun -machinefile hostfile -n ${nmpinodes[i]} -ppn ${npernode[i]} ${exec_script}" >> pbs/${exec_prefix}.pbs
		qsub -N ${exec_prefix} -l nodes=${nnodes[i]}:ppn=${ncores[i]},walltime=6:30:00 pbs/${exec_prefix}.pbs
	done
done

for i in {0..4}
do
	for j in {2..3}
	do
		# executing script
		exec_prefix="suc${nscenarios[i]}_${solver[j]}${cuts[j]}_oneiter"
		exec_script="julia suc_run_oneiter.jl ${nscenarios[i]} ${solver[j]} ${cuts[j]} ${exec_prefix}"

		cp pbs/template.pbs pbs/${exec_prefix}.pbs
		echo "mpirun -machinefile hostfile -n ${nmpinodes[i]} -ppn ${npernode[i]} ${exec_script}" >> pbs/${exec_prefix}.pbs
		qsub -N ${exec_prefix} -l nodes=${nnodes[i]}:ppn=${ncores[i]},walltime=6:30:00 pbs/${exec_prefix}.pbs
	done
done
