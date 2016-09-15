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
	"sslp_10_50_500"
	"sslp_10_50_1000"
	"sslp_10_50_2000"
)


for p in ${problems[@]}
do
	# executing script
	exec_prefix="${p}_bd"
	exec_args="${SIPLIB_DIR}/${p}"
	exec_script="julia smps_bd.jl ${exec_args} ${exec_prefix}"
	# execute!
	cp pbs/template.pbs pbs/${exec_prefix}.pbs
	echo "${exec_script}" >> pbs/${exec_prefix}.pbs
	qsub -N ${exec_prefix} -l nodes=1:ppn=1,walltime=1:00:00 pbs/${exec_prefix}.pbs
done
