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

        cp pbs/template.pbs pbs/${exec_prefix}.pbs
        echo "mpirun -machinefile hostfile -n ${nmpinodes[i]} -ppn ${npernode[i]} ${exec_script}" >> pbs/${exec_prefix}.pbs
        qsub -N ${exec_prefix} -l nodes=${nnodes[i]}:ppn=${ncores[i]},walltime=6:30:00 pbs/${exec_prefix}.pbs
done

for i in {0..3}
do
        # executing script
        exec_prefix="suctime5max${nperiods[i]}"
        exec_script="julia suc_run_time_dd.jl ${nperiods[i]} 1"

        cp pbs/template.pbs pbs/${exec_prefix}.pbs
        echo "mpirun -machinefile hostfile -n ${nmpinodes[i]} -ppn ${npernode[i]} ${exec_script}" >> pbs/${exec_prefix}.pbs
        qsub -N ${exec_prefix} -l nodes=${nnodes[i]}:ppn=${ncores[i]},walltime=6:30:00 pbs/${exec_prefix}.pbs
done
