#!/bin/bash

if [ -z "$DSP_INC" ]; then
	echo "Need to set DSP_INC; please read INSTALL file"
	exit 1
fi

SIPLIB_DIR="${DSP_INC}/../examples/smps"

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

for p in ${problems[@]}
do
	for i in {0..2}
	do
		# executing script
		exec_prefix="${p}_s${solvers[i]}"
		exec_args="${SIPLIB_DIR}/${p} ${solvers[i]} ${fcut[i]} ${ocut[i]} ${ddlog} ${iterlim} ${wt1}"
		exec_script="julia smps_dd.jl ${exec_args} ${exec_prefix}"
		mpirun -n ${nmpinodes} ${exec_script}
	done
done
