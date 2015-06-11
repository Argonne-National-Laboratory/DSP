#!/bin/bash

if [ -z "$DSP_INC" ]; then
	echo "Need to set DSP_INC; please read INSTALL file"
	exit 1
fi

SIPLIB_DIR="${DSP_INC}/../examples/smps"

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
	"sslp_10_50_500"
	"sslp_10_50_1000"
	"sslp_10_50_2000"
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
