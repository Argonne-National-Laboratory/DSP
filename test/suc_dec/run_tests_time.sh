#!/bin/bash

# algorithmic parameters
nperiods=(12 24 48 96)

# mpirun setting
nmpinodes=(2 4 8 16)
npernode=(2 4 8 16)

# cluster setting
nnodes=(1 1 1 1)
ncores=(2 4 8 16)

for i in {0..3}
do
        # executing script
        exec_prefix="suctime${nperiods[i]}"
        exec_script="julia suc_run_time_dd.jl ${nperiods[i]} 0"

        mpirun -n ${nmpinodes[i]} ${exec_script} > output/${exec_prefix}.out
done

for i in {0..3}
do
        # executing script
        exec_prefix="suctime5max${nperiods[i]}"
        exec_script="julia suc_run_time_dd.jl ${nperiods[i]} 1"

        mpirun -n ${nmpinodes[i]} ${exec_script} > output/${exec_prefix}.out
done
