#!/bin/bash

# algorithmic parameters
nscenarios=(4 8 16 32 64)
solver=2
cuts=1

# mpirun setting
nmpinodes=(4 8 16 32 64)
npernode=(4 8 16 16 16)

for i in {0..4}
do
        # executing script
        exec_prefix="suc${nscenarios[i]}-${solver}${cuts}"
        exec_script="julia suc_run.jl ${nscenarios[i]} ${solver} ${cuts} ${exec_prefix}"

        mpirun -n ${nmpinodes[i]} ${exec_script}
done
