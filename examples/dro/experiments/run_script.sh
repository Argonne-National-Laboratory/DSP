#!/bin/bash

#SBATCH --partition=bdwall
#SBATCH --account=NEXTGENOPT
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=32
#SBATCH --time=4:00:00
#SBATCH --mail-user=kimk@anl.gov
#SBATCH --mail-type=END

export I_MPI_FABRICS=shm:tmi
export OMP_NUM_THREADS=1

DSP=../../../build/bin/runDsp
K=$1
NP=64

echo "##########"
echo "# K=${K} #"
echo "##########"
echo ""

srun -n $NP $DSP --algo drdd --wassnorm 2.0 --wasseps 100.0 --smps drdcap_332_10_${K} --param params.txt > outputs/drdcap_332_${K}_100_${NP}.txt

# Run the following in the command line
# sbatch -J drdcap_scale_500 -o drdcap_scale_$(date +%Y%m%d%H%M%S).out run_script_scaling.sh 500

# Useful
# squeue -u kimk
# srun --pty -p bdwall -t 0:10:00 /bin/zsh
# sacct --format="JobID,JobName%20,NNodes,NTasks,CPUTime,MaxRSS,MaxRSSNode,MaxRSSTask,ExitCode" --units=G
# sstat --format=AveCPU,AvePages,MaxRSS,AveRSS,AveVMSize,JobID -j 1971639 --allsteps