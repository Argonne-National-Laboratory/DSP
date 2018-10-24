#!/bin/bash

maxcore=36
instances=( 
	"dcap233_200" "dcap233_300" "dcap233_500" 
	"dcap243_200" "dcap243_300" "dcap243_500" 
	"dcap332_200" "dcap332_300" "dcap332_500" 
	"dcap342_200" "dcap342_300" "dcap342_500" 
	"sslp_5_25_50" "sslp_5_25_100"
	"sslp_10_50_50" "sslp_10_50_100" "sslp_10_50_500" "sslp_10_50_1000" "sslp_10_50_2000"
	)
methods=( "BNP" "CS" "CSDW" )
heur=( "without" "with" )

small_instances=( "sslp_15_45_5" "sslp_15_45_10" "sslp_15_45_15" )
nprocs=( "5" "10" "15" )

for w in "${heur[@]}"
do
	for m in "${methods[@]}"
	do
		for i in "${instances[@]}"
		do
			echo "mpiexec -np $maxcore ../build/bin/runDsp --smps ../examples/smps/$i --param params/params-$m-$w-rounding.txt --soln solns/$i.$m.$w.sol.txt"
		done

		# Note: The number of processes should be no more than the number of scenarios
		for i in ${!small_instances[*]}
		do
			intance_name=${small_instances[$i]}
			echo "mpiexec -np ${nprocs[$i]} ../build/bin/runDsp --smps ../examples/smps/$intance_name --param params/params-$m-$w-rounding.txt --soln solns/$intance_name.$m.$w.sol.txt"
		done
	done
done
