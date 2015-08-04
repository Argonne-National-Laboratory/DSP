#!/bin/bash

# algorithmic parameters
nclusters=(2 3 4)

# mpirun setting
nmpinodes=(2 3 4)
npernode=(2 3 4)

# cluster setting
nnodes=(1 1 1)
ncores=(2 3 4)

for i in {0..2}
do
        # executing script
        exec_prefix="sucnetwork${nclusters[i]}"
        exec_script="julia suc_run_network.jl DD ${nclusters[i]}"

        mpirun -n ${nmpinodes[i]} ${exec_script} > output/${exec_prefix}.out
done
