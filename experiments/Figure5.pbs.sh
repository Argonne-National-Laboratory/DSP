#!/bin/bash

SIPLIB_DIR="./smps"

problems=(
	"sslp_15_45_15"
)

# algorithmic parameters
solvers=(3 0 2)
fcut=(-1 -1 1)
ocut=(-1 -1 1)
ddlog=1
iterlim=100
wt1=21600

# mpirun setting
nmpinodes=16
npernode=16

# cluster setting
nnodes=1
ncores=16

for p in ${problems[@]}
do
	for i in {0..2}
	do
		# executing script
		exec_prefix="${p}_s${solvers[i]}"
		exec_args="${SIPLIB_DIR}/${p} ${solvers[i]} ${fcut[i]} ${ocut[i]} ${ddlog} ${iterlim} ${wt1}"
		exec_script="julia smps_dd.jl ${exec_args} ${exec_prefix}"
		cp pbs/template.pbs pbs/${exec_prefix}.pbs
		echo "mpirun -machinefile hostfile -n ${nmpinodes} -ppn ${npernode} ${exec_script}" >> pbs/${exec_prefix}.pbs
		qsub -N ${exec_prefix} -l nodes=${nnodes}:ppn=${ncores},walltime=6:30:00 pbs/${exec_prefix}.pbs

	done
done
