#!/bin/bash

# algorithmic parameters
nscenarios=(4 8 16 32 64)

for i in {0..4}
do
	# executing script
	exec_prefix="suc${nscenarios[i]}-de"
	exec_script="julia suc_run_de.jl ${nscenarios[i]} ${exec_prefix}"

	cp pbs/template.pbs pbs/${exec_prefix}.pbs
	echo "mpirun -n 1 ${exec_script}" >> pbs/${exec_prefix}.pbs
	qsub -N ${exec_prefix} -l nodes=1:ppn=1,walltime=6:30:00 pbs/${exec_prefix}.pbs
done
