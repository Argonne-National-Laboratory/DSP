#!/bin/bash

SIPLIB_DIR="./smps"

problems=(
	"dcap233_200"
	"dcap233_300"
	"dcap233_500"
	"dcap243_200"
	"dcap243_300"
	"dcap243_500"
	"dcap332_200"
	"dcap332_300"
	"dcap332_500"
	"dcap342_200"
	"dcap342_300"
	"dcap342_500"
	"sslp_5_25_50"
	"sslp_5_25_100"
	"sslp_15_45_5"
	"sslp_15_45_10"
	"sslp_15_45_15"
	"sslp_10_50_50"
	"sslp_10_50_100"
)

# algorithmic parameters
solvers=(3 0 2)
fcut=(-1 -1 1)
ocut=(-1 -1 1)
ddlog=0
iterlim=10000
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
		# execute!
		cp pbs/template.pbs pbs/${exec_prefix}.pbs
		echo "mpirun -machinefile hostfile -n ${nmpinodes} -ppn ${npernode} ${exec_script}" >> pbs/${exec_prefix}.pbs
		qsub -N ${exec_prefix} -l nodes=${nnodes}:ppn=${ncores},walltime=6:30:00 pbs/${exec_prefix}.pbs
	done
done
