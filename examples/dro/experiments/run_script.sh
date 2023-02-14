#!/bin/bash

#SBATCH --partition=bdwall
#SBATCH --account=NEXTGENOPT
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=20
#SBATCH --time=2:00:00
#SBATCH --mail-user=kimk@anl.gov
#SBATCH --mail-type=END

export I_MPI_FABRICS=shm:tmi
export OMP_NUM_THREADS=1

DSP=../../../build/bin/runDsp
INSTANCES=( "233" "243" "332" "342" )
K=( "20" "50" "100" "200" "300" )
EPS=( "1" "100" "500" "1000" )

for i in "${INSTANCES[@]}"
do
    for k in "${K[@]}"
    do
        for e in "${EPS[@]}"
        do
            $DSP --algo drdd --wassnorm 2.0 --wasseps ${e}.0 --smps drdcap_${i}_10_${k} --param params.txt > outputs/drdcap_${i}_${k}_${e}.txt &
        done
    done
done

wait

# Run the following in the command line
# sbatch -J drdcap_500 -o drdcap_$(date +%Y%m%d%H%M%S).out run_script.sh 500

# Useful
# squeue -u kimk
# srun --pty -p bdwall -t 0:10:00 /bin/zsh
# sacct --format="JobID,JobName%20,NNodes,NTasks,CPUTime,MaxRSS,MaxRSSNode,MaxRSSTask,ExitCode" --units=G
# sstat --format=AveCPU,AvePages,MaxRSS,AveRSS,AveVMSize,JobID -j 1971639 --allsteps