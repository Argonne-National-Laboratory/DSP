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

        cp pbs/template.pbs pbs/${exec_prefix}.pbs
        echo "mpirun -machinefile hostfile -n ${nmpinodes[i]} -ppn ${npernode[i]} ${exec_script}" >> pbs/${exec_prefix}.pbs
        qsub -N ${exec_prefix} -l nodes=${nnodes[i]}:ppn=${ncores[i]},walltime=6:30:00 pbs/${exec_prefix}.pbs
done
