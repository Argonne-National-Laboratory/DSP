#!/bin/bash

SIPLIB_DIR="./smps"
siplib="sslp_10_50_500"

cut=1
maxeval=100

# mpirun setting
nmpinodes=(5 10 25 50 100 250)
npernode=(5 10 16 16 16 16)
wtimes=(1800 1800 1800 1800 3600 3600)
qsub_wtimes=("0:30:00" "0:30:00" "0:30:00" "0:30:00" "1:00:00" "1:00:00")

# cluster setting
nnodes=(1 1 2 4 7 16)
ncores=(5 10 16 16 16 16)
	
for i in {4..5}
do
	# executing script
	exec_prefix="${siplib}_scale${nmpinodes[i]}"
	exec_args="${SIPLIB_DIR}/${siplib} ${wtimes[i]} ${cut} ${maxeval}"
	exec_script="julia smps_dd_scale.jl ${exec_args} ${exec_prefix}"
	# execute!
	cp pbs/template.pbs pbs/${exec_prefix}.pbs
	echo "mpirun -machinefile hostfile -n ${nmpinodes[i]} -ppn ${npernode[i]} ${exec_script}" >> pbs/${exec_prefix}.pbs
	qsub -N ${exec_prefix} -l nodes=${nnodes[i]}:ppn=${ncores[i]},walltime=${qsub_wtimes[i]} pbs/${exec_prefix}.pbs
done
